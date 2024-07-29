/*
 * Multiply a, b in (Z/p)[x]
 * We require that:
 *  degree(a) + degree(b) < n; and
 *  the arrays a and b both have length at least n.
 * 
 * Store the result in a, but return false if no primitive n-th root of unity
 * could be found in Z/p.
 */
bool fast_multiply(int *const a, int *const b, int n, int p);
