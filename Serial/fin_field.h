#ifdef __cplusplus
    extern "C" {
#endif

#define PRIME_LOG_2  16


// Compute a*b mod p
int mulmod(int a, int b);

// Compute a + b mod p
int addmod(int a, int b);

// Compute a - b mod p
int submod(int a, int b);

// compute a^n mod p
int powmod(int a, int n);

// Generate a random element of Z/p
int rand_elt(void);

// Let p be a Fermat prime. Return a primitive nth root of unity in Z/p if one
// exists. If not, return 0.
int fermat_primitive_root(int n);

// Find the (multiplicative) inverse of n in Z/p
int invmod(int n);

#ifdef __cplusplus
    }
#endif
