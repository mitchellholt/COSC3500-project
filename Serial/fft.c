#include "fft.h"
#include "fin_field.h"
#include <stdint.h>
#include <stdlib.h>


// EPIC implementation that uses n/2 field multiplications and O(n) array
// operations. TODO - can we make this n accesses?
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
    // TODO can we force w to stay in L1?
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


int *fft1(int *const coefficients, int n, const int *const w, int p) {
    return NULL;
}


int *fft2(int *const coefficients, int n, const int *const w, int p) {
    return NULL;
}
