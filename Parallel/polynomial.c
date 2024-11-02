#include <mm_malloc.h>
#include <stdlib.h>
#include <string.h>

#include "polynomial.h"
#include "fft.h"
#include "fin_field.h"

#define ALIGN           64
#define CACHE_LINE_INTS 4


void fast_multiply(int *const a, int *const b, int omega, int n, int p) {
    int *w = _mm_malloc(sizeof(int) * n, ALIGN);
    primitive_root_powers(w, n, omega, p);
    fft2(a, n, w, p, 0);
    fft2(b, n, w, p, 0);

    const int vec_n = n/CACHE_LINE_INTS;
    __m128i *a_vec = (__m128i*)a;
    __m128i *b_vec = (__m128i*)b;
    __m128i inv = _mm_set1_epi32(invmod(n, p));

    for (int i = 0; i < vec_n; i++) { // NOTE: n is always divisible by 8
        __m128i b_block = b_vec[i];
        a_vec[i] = vec_mulmod(vec_mulmod(a_vec[i], inv, p), b_block, p);
    }

    omega = invmod(omega, p);
    primitive_root_powers(w, n, omega, p);
    fft1(a, n, w, p, 0);

    _mm_free(w);
}


void rand_poly(int *buffer, int degree, int p) {
    for (int i = 0; i < degree; i++) {
        buffer[i] = rand_elt(p);
    }
}


int degree(int n, int *poly) {
    int deg = n - 1;
    for(; deg > 0 && !poly[deg]; deg--);
    return deg;
}
