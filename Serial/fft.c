#include <stdint.h>
#include <mm_malloc.h>
#include "fft.h"
#include "fin_field.h"


// New algrithm does twice as many multiplications, but buffer should be in
// cache the whole time.
void primitive_root_powers(int *buffer, int n, int omega) {
    if (n == 0) return;

    int n2 = n/2;
    int om = 1;
    for (int i = 0; i < n2; i++) {
        buffer[i] = om;
        om = mulmod(om, omega);
    }

    primitive_root_powers(buffer + n, n2, mulmod(omega, omega));
}


void fft1(int *const coefficients, int n, const int *const w) {
    if (n == 1) return;
    const int n2 = n/2;
    fft1(coefficients, n2, w + n2);
    fft1(coefficients + n2, n2, w + n2);
    int s,t;
    for (int i = 0; i < n2; i++) {
        s = coefficients[i];
        t = mulmod(w[i], coefficients[n2 + i]);
        coefficients[i] = addmod(s, t);
        coefficients[n2 + i] = submod(s, t);
    }
    return;
}


void fft2(int *const coefficients, int n, const int *const w) {
    if (n == 1) return;
    const int n2 = n/2;
    int s,t;
    for (int i = 0; i < n2; i++) {
        s = addmod(coefficients[i], coefficients[n2 + i]);
        t = submod(coefficients[i], coefficients[n2 + i]);
        coefficients[i] = s;
        coefficients[n2 + i] = mulmod(t, w[i]);
    }
    fft2(coefficients, n2, w + n2);
    fft2(coefficients + n2, n2, w + n2);
    return;
}
