#include <stdint.h>
#include <mm_malloc.h>
#include "fft.h"
#include "fin_field.h"


// Note that this is fast IF we can fit all of w into cache.
// What happens otherwise?
void primitive_root_powers(int *buffer, int n, int omega) {
    // fill in first half of the array (multiplications should ONLY be done here)
    int k = n/2; // number of values to write here
    int om = omega;
    buffer[0] = 1;
    for (int i = 1; i < k; i++) {
        buffer[i] = om;
        om = mulmod(om, omega);
    }

    // fill in the rest of the values, using those already in the start of the array
    // TODO: use epic cache stuff to try to keep the parts of w we need in cache
    uint16_t len = k; len >>= 1; // number of values to write in the current iteration
    uint16_t shift = 1; // we write the (multiple)th powers of omega in the current iteration
    // From here, k is the start index of the current block
    while (len) {
        for (uint16_t i = 0; i < len; i++) {
            buffer[k + i] = buffer[i << shift];
        }
        k += len;
        len >>= 1;
        shift++;
    }
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
