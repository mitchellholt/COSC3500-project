#include "fft.h"
#include "fin_field.h"
#include <stdint.h>
#include <stdlib.h>


// Note that this is fast IF we can fit all of w into cache.
// What happens otherwise?
int *primitive_root_powers(int n, int omega, int p) {
    int *w = malloc(sizeof(int) * n); // no need for null-termination

    // fill in first half of the array (multiplications should ONLY be done here)
    int k = n/2; // number of values to write here
    int om = omega;
    w[0] = 1;
    for (int i = 1; i < k; i++) {
        w[i] = om;
        om = mulmod(om, omega, p);
    }

    // fill in the rest of the values, using those already in the start of the array
    // TODO: use epic cache stuff to try to keep the parts of w we need in cache
    uint16_t len = k; len >>= 1; // number of values to write in the current iteration
    uint16_t shift = 1; // we write the (multiple)th powers of omega in the current iteration
    // From here, k is the start index of the current block
    while (len) {
        for (uint16_t i = 0; i < len; i++) {
            w[k + i] = w[i << shift];
        }
        k += len;
        len >>= 1;
        shift++;
    }

    return w;
}


void fft1(int *const coefficients, int n, const int *const w, int p) {
    if (n == 1) return;
    const int n2 = n/2;
    fft1(coefficients, n2, w + n2, p);
    fft1(coefficients + n2, n2, w + n2, p);
    int s,t;
    for (int i = 0; i < n2; i++) {
        s = coefficients[i];
        t = mulmod(w[i], coefficients[n2 + i], p);
        coefficients[i] = addmod(s, t, p);
        coefficients[n2 + i] = submod(s, t, p);
    }
    return;
}


void fft2(int *const coefficients, int n, const int *const w, int p) {
    if (n == 1) return;
    const int n2 = n/2;
    int s,t;
    for (int i = 0; i < n2; i++) {
        s = addmod(coefficients[i], coefficients[n2 + i], p);
        t = submod(coefficients[i], coefficients[n2 + i], p);
        coefficients[i] = s;
        coefficients[n2 + i] = mulmod(t, w[i], p);
    }
    fft2(coefficients, n2, w + n2, p);
    fft2(coefficients + n2, n2, w + n2, p);
    return;
}
