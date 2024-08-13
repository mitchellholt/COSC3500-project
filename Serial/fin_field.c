#include "fin_field.h"
#include <stdint.h>
#include <stdlib.h>


inline int mulmod(int a, int b) {
    // Since the prime is 2^16 + 1, we can compute q,r as solutions to the
    // Euclidean division
    //  a*b = q*2^16 + r
    // Then a*b mod (2^16 + 1) = r - q
    uint64_t t = ((uint64_t)a) * ((uint64_t)b); // a <= 2^16, b <= 2^16, so ab <= 2^32. Need a 64 bit integer
    int q = (int)(t >> PRIME_LOG_2);
    int r = (int)(t & 0xFFFF); // r is the least significant 16 bits
    return submod(r, q);
}


inline int addmod(int a, int b) {
    int t = a - FERMAT_PRIME + b;
    t += (t >> 31) & FERMAT_PRIME;
    return t;
}


inline int submod(int a, int b) {
    int t = a - b;
    t += (t >> 31) & FERMAT_PRIME;
    return t;
}


inline int powmod(int a, int n) {
    int t = 1;
    for (; n > 0; n--) t = mulmod(t, a);
    return t;
}


inline int rand_elt(void) {
    return rand() % FERMAT_PRIME;
}


// Assume p is a Fermat prime and n is a power of 2.
// Moreover, assume that p <= 2^16 + 1
int fermat_primitive_root(int n) {
    if (n >= FERMAT_PRIME) return 0;

    // 2 is the ONLY prime divisor of FERMAT_PRIME - 1
    int candidate = 0; // candidate (FERMAT_PRIME - 1)th root of unity
    while (!candidate) { // expected to run through this loop 3 times
        candidate = (rand() % (FERMAT_PRIME - 2)) + 2;
        if (powmod(candidate, (FERMAT_PRIME - 1)/2) == 1) {
            candidate = 0;
            continue;
        }
    }

    return powmod(candidate, (FERMAT_PRIME - 1)/n);
}


// PERF: wow this kind of sucks. We could use repeated squaring to speed up
// powmod or implement XGCD
inline int invmod(int n) {
    return powmod(n, FERMAT_PRIME - 2);
}
