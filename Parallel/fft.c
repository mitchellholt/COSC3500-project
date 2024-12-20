#include <stdint.h>
#include <mm_malloc.h>
#include <string.h>
#include <xmmintrin.h>
#include <omp.h>

#include "fft.h"
#include "fin_field.h"

#define CACHE_LINE_INTS    4
#define PARALLEL_MAX_DEPTH 2
#define PARALLEL_MIN_SIZE  16


// Note that this is fast IF we can fit all of w into cache.
// What happens otherwise?
void primitive_root_powers(int *buffer, int n, int omega, int p) {
    // fill in first half of the array (multiplications should ONLY be done here)
    int k = n/2; // number of values to write here
    int om = omega;
    buffer[0] = 1;
    for (int i = 1; i < k; i++) {
        buffer[i] = om;
        om = mulmod(om, omega, p);
    }

    // fill in the rest of the values, using those already in the start of the array
    uint16_t len = k; len >>= 1; // number of values to write in the current iteration
    uint16_t shift = 1; // we write the (multiple)th powers of omega in the current iteration
    // From here, k is the start index of the current block
    // PERF: This is an actual cache miss machine
    while (len) {
        for (uint16_t i = 0; i < len; i++) {
            buffer[k + i] = buffer[i << shift];
        }
        k += len;
        len >>= 1;
        shift++;
    }
}


static inline __m128i fft1_base_case(__m128i coeffs, __m128i w, int p) {
    // compute s,t vectors
    __m128i t = _mm_shuffle_epi32(coeffs, _MM_SHUFFLE(3, 2, 3, 2));
    __m128i s = _mm_shuffle_epi32(coeffs, _MM_SHUFFLE(1, 0, 1, 0));
    t = vec_mulmod(t, w, p);

    // store the result
    __m128i adds = vec_addmod(s, t, p);
    __m128i subs = vec_submod(s, t, p);
    return _mm_blend_epi32(adds, subs, 0xC); // 0b0011
}


// Based on Michael Monagan's code
// https://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA21/FFTnoperm.pdf
// WARN: Only works when n >= 4
void fft1(int *const coefficients, int n, const int *const w, int p, int depth) {
    __m128i *coeffs_vec = (__m128i*)coefficients;
    __m128i *w_vec = (__m128i*)w;

    if (n == 4) { // write this explicitly so that 4|n later
        __m128i w_sparse = _mm_shuffle_epi32(*w_vec, _MM_SHUFFLE(2, 0, 2, 0));
        __m128i w_lower = _mm_shuffle_epi32(*w_vec, _MM_SHUFFLE(1, 0, 1, 0));
        __m128i coeffs = _mm_shuffle_epi32(*coeffs_vec, _MM_SHUFFLE(3, 1, 2, 0));

        coeffs = _mm_shuffle_epi32(
                fft1_base_case(coeffs, w_sparse, p),
                _MM_SHUFFLE(3, 1, 2, 0));
        *coeffs_vec = fft1_base_case(coeffs, w_lower, p);
        return;
    }

    const int vec_n2 = n/(2*CACHE_LINE_INTS);
    const int n2 = CACHE_LINE_INTS*vec_n2;
    fft1(coefficients,      n2, w + n2, p, depth + 1);
    fft1(coefficients + n2, n2, w + n2, p, depth + 1);

    for (int i = 0; i < vec_n2; i++) {
        __m128i b = coeffs_vec[vec_n2 + i];
        __m128i t = vec_mulmod(w_vec[i], b, p);
        __m128i s = coeffs_vec[i];

        coeffs_vec[i] = vec_addmod(s, t, p);
        coeffs_vec[vec_n2 + i] = vec_submod(s, t, p);
    }
    return;
}


static inline __m128i fft2_base_case(__m128i coeffs, __m128i w, int p) {
    // things to subtract/add
    __m128i ops = _mm_shuffle_epi32(coeffs, _MM_SHUFFLE(1, 0, 1, 0));
    __m128i opa = _mm_shuffle_epi32(coeffs, _MM_SHUFFLE(3, 2, 3, 2));

    __m128i t = vec_submod(coeffs, ops, p);
    __m128i s = vec_addmod(coeffs, opa, p);
    t = vec_mulmod(t, w, p);

    return _mm_blend_epi32(s, t, 0xC);
}


// Based on Michael Monagan's code 
// https://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA21/FFTnoperm.pdf
void fft2(int *const coefficients, int n, const int *const w, int p, int depth) {
    __m128i *coeff_vec = (__m128i*)coefficients;
    __m128i *w_vec = (__m128i*)w;

    if (n == 4) { // BASE CASE - n = 2, n = 1 have been manually inlined.
        __m128i w2 = _mm_shuffle_epi32(*w_vec, _MM_SHUFFLE(2, 2, 2, 2));
        __m128i w_lower = _mm_shuffle_epi32(*w_vec, _MM_SHUFFLE(1, 0, 1, 0));
        __m128i coeffs = *coeff_vec;

        coeffs = _mm_shuffle_epi32(
                fft2_base_case(coeffs, w_lower, p),
                _MM_SHUFFLE(3, 1, 2, 0));

        *coeff_vec = _mm_shuffle_epi32(
                fft2_base_case(coeffs, w2, p),
                _MM_SHUFFLE(3, 1, 2, 0));

        return;
    }

    // NOTE: coefficients and w should always be aligned along a 4*4=16 byte
    // boundary
    const int vec_n2 = n/(2*CACHE_LINE_INTS);

    for (int i = 0; i < vec_n2; i++) {
        __m128i b = coeff_vec[vec_n2 + i];
        __m128i t = vec_submod(coeff_vec[i], b, p);

        coeff_vec[i] = vec_addmod(coeff_vec[i], b, p);
        coeff_vec[vec_n2 + i] = vec_mulmod(t, w_vec[i], p);
    }

    const int n2 = CACHE_LINE_INTS*vec_n2;

    if (depth <= PARALLEL_MAX_DEPTH) {
        #pragma omp parallel sections
        {
            #pragma omp section
            { fft2(coefficients,      n2, w + n2, p, depth + 1); }
            #pragma omp section
            { fft2(coefficients + n2, n2, w + n2, p, depth + 1); }
        }
    } else {
        fft2(coefficients,      n2, w + n2, p, depth + 1);
        fft2(coefficients + n2, n2, w + n2, p, depth + 1);
    }

    return;
}
