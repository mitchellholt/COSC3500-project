#ifdef __cplusplus
    extern "C" {
#endif
/*
 * Multiply a, b in (Z/p)[x]
 * We require that:
 *  degree(a) + degree(b) < n,
 *  n is a power of 2,
 *  the arrays a and b both have length at least n, and
 *  omega is a primitive nth root of unity in Z/p
 * 
 * Store the result in a, but return false if no primitive n-th root of unity
 * could be found in Z/p.
 */
void fast_multiply(int *const a, int *const b, int omega, int n, int p);


/*
 * Create a random polynomial in (Z/p)[x] with at most the specified degree and
 * store it in buffer
 *
 * The caller must ensure that buffer is large enough to hold degree ints.
 */
void rand_poly(int *buffer, int degree, int p);

/*
 * Return the degree of the input polynomial, which is stored in a buffer with n
 * elements
 */
int degree(int n, int *poly);

#ifdef __cplusplus
    }
#endif
