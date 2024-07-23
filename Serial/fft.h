/* Compute an array of length n of powers of omega modulo p. We require that
 * omega is a primitive nth root of unity.
 *
 * The format of this array is identical
 * to that in Michael Monagan's implementnation
 * (https://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA21/FFTnoperm.pdf),
 * with the exception that the final element of the array is left uninitialised
 */
int *primitive_root_powers(int n, int omega, int p);


/* In-place FFT Algorithm 1 (Geddes, Czapor, Labahn section 4.7)
 * We require that:
 *  - n is a power of two,
 *  - p is prime,
 *  - coefficients has length n,
 *  - w is a precomputed array of powers of some primitive n-th root of unity
 *    (computed by primitive_root_powers) in Z/p.
 *
 * The result is stored in coefficients; no new memory is allocated. A pointer
 * to the coefficient array is also returned.
 */
int *fft1(int *coefficients, int n, const int *w, int p);


/* In-place FFT Algorithm 2 (Gerhard and von zur Gathen section 8.2)
 * We require that:
 *  - n is a power of two,
 *  - p is prime,
 *  - coefficients has length n,
 *  - w is a precomputed array of powers of some primitive n-th root of unity
 *    (computed by primitive_root_powers) in Z/p.
 *
 * The result is stored in coefficients; no new memory is allocated. A pointer
 * to the coefficient array is also returned.
 */
int *fft2(int *coefficients, int n, const int *w, int p);
