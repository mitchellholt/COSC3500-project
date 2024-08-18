#include <mm_malloc.h>
#include <stdbool.h>
#include <stdlib.h>

#include "polynomial.h"
#include "fft.h"
#include "fin_field.h"

#define ALIGN 64


// PERF: can we somehow cache w(omega) and w(omega^{-1})?
void fast_multiply(int *const a, int *const b, int omega, int n, int p) {
    int *w = _mm_malloc(sizeof(int) * n, ALIGN);
    primitive_root_powers(w, n, omega, p);
    fft2(a, n, w, p);
    fft2(b, n, w, p);

    // PERF: Consider cache blocks
    // TODO: this step is begging for SIMD
    int n_inv = invmod(n, p);
    for (int i = 0; i < n; i++) {
        // NOTE: We can do the multiplication here (instead of after fft1)
        // because FFT is linear. Probably improves CACHE locality
        a[i] = mulmod(n_inv, mulmod(a[i], b[i], p), p);
    }

    omega = invmod(omega, p);
    // PERF: surely there is an efficient way to do this WITHOUT needing to create a new array?
    primitive_root_powers(w, n, omega, p);
    fft1(a, n, w, p);

    free(w);
}


void rand_poly(int *buffer, int degree, int p) {
    for (int i = 0; i < degree; i++) {
        buffer[i] = rand_elt(p);
    }
}
