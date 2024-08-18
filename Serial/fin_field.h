#ifdef __cplusplus
    extern "C" {
#endif

// Compute a*b mod p
int mulmod(int a, int b, int p);

// Compute a + b mod p
int addmod(int a, int b, int p);

// Compute a - b mod p
int submod(int a, int b, int p);

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
