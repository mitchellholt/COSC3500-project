#include "fin_field.h"
#include <stdlib.h>

#define LONG long long int // no need to expose this


inline int mulmod(int a, int b, int p) {
    int t = (LONG) a * b % p;
    return t;
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


inline int powmod(int a, int n, int p) {
    int t = 1;
    for (; n > 0; n--) t = mulmod(t, a, p);
    return t;
}


int primitive_root(int n, int p) {
    // Section 4.8 Of Algorithms for Computer Algebra (Geddes, Czapor, Labahn)
    if ((p - 1) % n) return 0; // Theorem 4.3

    // for now, set prime_divisors to be all primes less than p - 1
    int *prime_divisors = malloc(sizeof(int) * (p - 1));
    int num_prime_divisors = 0;
    for (int k = 2; k < p - 1; k++) {
        int j;
        for (j = 0; j < num_prime_divisors; j++) {
            if (k % prime_divisors[j] == 0) {
                break;
            }
        }
        if (j == num_prime_divisors) {
            prime_divisors[num_prime_divisors++] = k;
        }
    }

    for (int i = 0; i < num_prime_divisors; i++) {
        if ((p - 1) % prime_divisors[i] != 0) prime_divisors[i] = 0;
    }

    int candidate = 0; // candidate (p - 1)th root of unity
    while (!candidate) { // expected to run through this loop 3 times
        candidate = (rand() % (p - 2)) + 2;
        for (int i = 0; i < num_prime_divisors; i++) {
            if (prime_divisors[i] && powmod(candidate, (p - 1)/prime_divisors[i], p) == 1) {
                candidate = 0;
                break;
            }
        }
    }

    free(prime_divisors);
    return powmod(candidate, (p - 1)/n, p);
}

// PERF: maybe this would be faster using XGCD?
inline int invmod(int n, int p) {
    return powmod(n, p - 2, p);
}
