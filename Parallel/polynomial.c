#include <mm_malloc.h>
#include <stdlib.h>
#include <string.h>

#include "polynomial.h"
#include "fft.h"
#include "fin_field.h"

#define ALIGN           64
#define CACHE_LINE_INTS 4

// int triple_mulmod(int a, int b, int c, int p);
// template<int N> void fn() {
//     POINTWISE_PRODUCT(N-1);
//     fn<N-1>();
// }
// template<> void fn<0>() {
// }


void fast_multiply(int *const a, int *const b, int omega, int n, int p) {
    int *w = _mm_malloc(sizeof(int) * n, ALIGN);
    primitive_root_powers(w, n, omega, p);
    fft2(a, n, w, p);
    fft2(b, n, w, p);

    // TODO: this step is begging for SIMD
    int b_block[CACHE_LINE_INTS];
    int n_inv = invmod(n, p);
    for (int i = 0; i < n; i += CACHE_LINE_INTS) { // NOTE: n is always divisible by 8
        memcpy(b_block, b + i, sizeof(int) * CACHE_LINE_INTS);
        a[i    ] = triple_mulmod(n_inv, a[i    ], b_block[0], p);
        a[i + 1] = triple_mulmod(n_inv, a[i + 1], b_block[1], p);
        a[i + 2] = triple_mulmod(n_inv, a[i + 2], b_block[2], p);
        a[i + 3] = triple_mulmod(n_inv, a[i + 3], b_block[3], p);
    }

    omega = invmod(omega, p);
    // PERF: surely there is an efficient way to do this WITHOUT needing to create a new array?
    primitive_root_powers(w, n, omega, p);
    fft1(a, n, w, p);

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
