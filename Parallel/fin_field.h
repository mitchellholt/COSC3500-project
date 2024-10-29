#ifdef __cplusplus
    extern "C" {
#endif

#include <immintrin.h>

// Compute a*b mod p
int mulmod(int a, int b, int p);

// Compute a*b*c mod p
int triple_mulmod(int a, int b, int c, int p);

// __m128i vec_triple_mulmod(__m128i x, __m128i y, __m128i z, int p);

// Compute the product of the packed 32-bit integer vectros x,y mod p.
// Require p = 2^16 + 1
__m128i vec_mulmod(__m128i x, __m128i y, int p);

// Compute a + b mod p
int addmod(int a, int b, int p);

// Pointwise addition (mod p) of 128 bit packed vectors of integers
__m128i vec_addmod(__m128i x, __m128i y, int p);

// Compute a - b mod p
int submod(int a, int b, int p);

__m128i vec_submod(__m128i x, __m128i y, int p);

// compute a^n mod p
int powmod(int a, int n, int p);

// Generate a random element of Z/p
int rand_elt(int p);

// Let p be a Fermat prime. Return a primitive nth root of unity in Z/p if one
// exists. If not, return 0.
int fermat_primitive_root(int n, int p);

// Return a primitive nth root of unity in Z/p if one exists. If not, return 0
// int primitive_root(int n, int p);

// Find the (multiplicative) inverse of n in Z/p
int invmod(int n, int p);

#ifdef __cplusplus
    }
#endif
