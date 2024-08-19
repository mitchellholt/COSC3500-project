#include "fin_field.h"
#include <stdint.h>
#include <stdlib.h>


inline int mulmod(int a, int b, int p) {
    int64_t t = ((int64_t)a * (int64_t)b) % p;
    return (int)t;
}


inline int triple_mulmod(int a, int b, int c, int p) {
    int64_t t = ((int64_t)a * (int64_t)b * (int64_t)c) % p;
    return (int)t;
}


inline int addmod(int a, int b, int p) {
    int t = a - p + b;
    t += (t >> 31) & p;
    return t;
}


inline int submod(int a, int b, int p) {
    int t = a - b;
    t += (t >> 31) & p;
    return t;
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
