#include "fin_field.h"

#include <stdint.h>
#include <stdlib.h>
#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>
#include <immintrin.h>


inline int mulmod(int a, int b, int p) {
    int64_t t = ((int64_t)a * (int64_t)b) % p;
    return (int)t;
}


// Input is packed 32-bit integers
inline __m128i vec_mulmod(__m128i x, __m128i y, int p) {
    // mask to get the 16 least significant bits
    __m128i mask = _mm_set1_epi64x(0xFFFF);
    // get the high (odd) values from the vectors
    __m128i high_x = _mm_shuffle_epi32(x, _MM_SHUFFLE(2, 3, 0, 1));
    __m128i high_y = _mm_shuffle_epi32(y, _MM_SHUFFLE(2, 3, 0, 1));
    // Compute the product in Z; ouput is 64 bit unsigned integers
    __m128i high_prod = _mm_mul_epu32(high_x, high_y);
    __m128i low_prod = _mm_mul_epu32(x, y);
    // Compute quotients and remainders
    __m128i high_quo = _mm_srli_epi64(high_prod, 16);
    __m128i low_quo = _mm_srli_epi64(low_prod, 16);
    __m128i high_rem = _mm_and_si128(high_prod, mask);
    __m128i low_rem = _mm_and_si128(low_prod, mask);
    // Pack into 32 bit integers
    high_quo = _mm_shuffle_epi32(high_quo, _MM_SHUFFLE(2, 3, 0, 1));
    high_rem = _mm_shuffle_epi32(high_rem, _MM_SHUFFLE(2, 3, 0, 1));
    __m128i quo = _mm_blend_epi32(high_quo, low_quo, 5); // 5 = 0b0101
    __m128i rem = _mm_blend_epi32(high_rem, low_rem, 5); // 5 = 0b0101
    // Finally, do the subtraction!
    return vec_submod(rem, quo, p);
}


inline int triple_mulmod(int a, int b, int c, int p) {
    int64_t t = ((int64_t)a * (int64_t)b * (int64_t)c) % p;
    return (int)t;
}


inline __m128i vec_triple_mulmod(__m128i x, __m128i y, __m128i z, int p) {
    __m128i zero = _mm_set1_epi64x(0);
    // mask to get the 16 least significant bits
    __m128i mask = _mm_set1_epi64x(0xFFFF);

    // get the high (odd) values from the vectors
    __m128i high_x = _mm_shuffle_epi32(x, _MM_SHUFFLE(2, 3, 0, 1));
    __m128i high_y = _mm_shuffle_epi32(y, _MM_SHUFFLE(2, 3, 0, 1));
    __m128i high_z = _mm_shuffle_epi32(z, _MM_SHUFFLE(2, 3, 0, 1));
    // Compute the product in Z; ouput is 64 bit unsigned integers
    __m128i high_prod = _mm_mul_epu32(high_x, high_y);
    __m128i low_prod = _mm_mul_epu32(x, y);
    __m128i z_uppers = _mm_blend_epi32(high_z, zero, 0xA); // 0b1010
    __m128i z_lowers = _mm_blend_epi32(z, zero, 0xA); // 0b1010
    high_prod = _mm_mullo_epi64(high_prod, z_uppers); // there won't be any overflow (but we still have to eat the cost)
    low_prod = _mm_mullo_epi64(low_prod, z_lowers);

    // Compute quotients and remainders
    __m128i high_quo = _mm_srli_epi64(high_prod, 16);
    __m128i low_quo = _mm_srli_epi64(low_prod, 16);
    __m128i high_rem = _mm_and_si128(high_prod, mask);
    __m128i low_rem = _mm_and_si128(low_prod, mask);
    // Pack into 32 bit integers
    high_quo = _mm_shuffle_epi32(high_quo, _MM_SHUFFLE(2, 3, 0, 1));
    high_rem = _mm_shuffle_epi32(high_rem, _MM_SHUFFLE(2, 3, 0, 1));
    __m128i quo = _mm_blend_epi32(high_quo, low_quo, 0x5); // 0b0101
    __m128i rem = _mm_blend_epi32(high_rem, low_rem, 0x5); // 0b0101
    // Finally, do the subtraction!
    return vec_submod(rem, quo, p);
}


inline int addmod(int a, int b, int p) {
    int t = a - p + b;
    t += (t >> 31) & p;
    return t;
}


// algorithm due to van der Hoeven and Lecerf (https://arxiv.org/pdf/1407.3383)
// PERF: see what happens if we up this bad boy to 256-bit vectors.
// May not be faster bc memory throughput (cringe)
inline __m128i vec_addmod(__m128i x, __m128i y, int p) {
    const __m128i p_vec = _mm_set1_epi32(p);
    __m128i a = _mm_sub_epi32(p_vec, y);
    __m128i b = _mm_cmpeq_epi32(_mm_max_epu32(x, a), x);
    __m128i c = _mm_andnot_si128(b, p_vec);
    return _mm_add_epi32(_mm_sub_epi32(x, a), c);
}


inline int submod(int a, int b, int p) {
    int t = a - b;
    t += (t >> 31) & p;
    return t;
}


// I wrote this myself
inline __m128i vec_submod(__m128i x, __m128i y, int p) {
    const __m128i p_vec = _mm_set1_epi32(p);
    __m128i a = _mm_sub_epi32(x, y);
    __m128i b = _mm_cmpeq_epi32(_mm_max_epu32(x, a), x); // overflow <==> this is 0
    __m128i c = _mm_andnot_si128(b, p_vec);
    return _mm_add_epi32(a, c);
}


// PERF: I reckon that this is almost always called with n = 2^k - 1, in which
// case we get worst case performance.
// TODO: Is there a way to improve this?
int powmod(int a, int n, int p) {
    if (n == 0) return 1;
    if (n == 1) return a;

    uint32_t log2n = 31 - __builtin_clz((uint32_t )n);
    int t = a;
    for (uint32_t i = 0; i < log2n; i++) {
        // loop invariant: t = a^(2^i)
        t = mulmod(t, t, p);
    }
    int tt = powmod(a, n - (int)(1 << log2n), p);
    return mulmod(t, tt, p);
}


inline int rand_elt(int p) {
    return rand() % p;
}


// Assume p is a Fermat prime and n is a power of 2.
// Moreover, assume that p <= 2^16 + 1
int fermat_primitive_root(int n, int p) {
    if (n >= p) return 0;

    // 2 is the ONLY prime divisor of p - 1
    int candidate = 0; // candidate (p - 1)th root of unity
    while (!candidate) { // expected to run through this loop 3 times
        candidate = (rand() % (p - 2)) + 2;
        if (powmod(candidate, (p - 1)/2, p) == 1) {
            candidate = 0;
            continue;
        }
    }

    return powmod(candidate, (p - 1)/n, p);
}


// PERF: wow this kind of sucks. We could use repeated squaring to speed up
// powmod or implement XGCD
inline int invmod(int n, int p) {
    return powmod(n, p - 2, p);
}
