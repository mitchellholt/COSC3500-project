#include <stdint.h>
#include <mm_malloc.h>
#include <string.h>
#include "fft.h"
#include "fin_field.h"

#define CACHE_LINE_INTS 4


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


// Based on Michael Monagan's code
// https://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA21/FFTnoperm.pdf
void fft1(int *const coefficients, int n, const int *const w, int p) {
    if (n == 1) return; // base case

    if (n == 2) { // write this explicitly so that 4|n later
        // recursive FFT calls with n=1 do nothing
        int s = coefficients[0];
        int t = mulmod(w[0], coefficients[1], p);
        coefficients[0] = addmod(s, t, p);
        coefficients[1] = submod(s, t, p);
        return;
    }

    if (n == 4) {
        fft1(coefficients, 2, w + 2, p);
        fft1(coefficients + 2, 2, w + 2, p);

        int s1 = coefficients[0];
        int s2 = coefficients[1];
        int t1 = mulmod(w[0], coefficients[2], p);
        int t2 = mulmod(w[1], coefficients[3], p);

        coefficients[0] = addmod(s1, t1, p);
        coefficients[1] = addmod(s2, t2, p);
        coefficients[2] = submod(s1, t1, p);
        coefficients[3] = submod(s2, t2, p);

        return;
    }

    const int n2 = n/2;
    fft1(coefficients, n2, w + n2, p);
    fft1(coefficients + n2, n2, w + n2, p);
    int s,t;

    int bbuf[CACHE_LINE_INTS];
    int wbuf[CACHE_LINE_INTS];
    for (int i = 0; i < n2; i += CACHE_LINE_INTS) { // NOTE: n2 is divisible by 4:
        memcpy(bbuf, coefficients + n2 + i, sizeof(int) * CACHE_LINE_INTS);
        memcpy(wbuf, w + i, sizeof(int) * CACHE_LINE_INTS);

        s = coefficients[i    ];
        t = mulmod(wbuf[0], bbuf[0], p);
        coefficients[i    ] = addmod(s, t, p);
        bbuf[0] = submod(s, t, p);

        s = coefficients[i + 1];
        t = mulmod(wbuf[1], bbuf[1], p);
        coefficients[i + 1] = addmod(s, t, p);
        bbuf[1] = submod(s, t, p);

        s = coefficients[i + 2];
        t = mulmod(wbuf[2], bbuf[2], p);
        coefficients[i + 2] = addmod(s, t, p);
        bbuf[2] = submod(s, t, p);

        s = coefficients[i + 3];
        t = mulmod(wbuf[3], bbuf[3], p);
        coefficients[i + 3] = addmod(s, t, p);
        bbuf[3] = submod(s, t, p);

        memcpy(coefficients + n2 + i, bbuf, sizeof(int) * CACHE_LINE_INTS);
    }

    return;
}


// Based on Michael Monagan's code 
// https://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA21/FFTnoperm.pdf
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
