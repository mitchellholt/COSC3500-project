#ifdef __cplusplus
    extern "C" {
#endif

/* Compute an array of length n of powers of omega modulo p. We require that
 * omega is a primitive nth root of unity. Store them in the buffer that is
 * passed to this function. The caller must ensure that the buffer is large
 * enough
 *
 * The format of this array is identical
 * to that in Michael Monagan's implementnation
 * (https://www.cecm.sfu.ca/~mmonagan/teaching/TopicsinCA21/FFTnoperm.pdf),
 * with the exception that the final element of the array is left uninitialised
 */
void primitive_root_powers(int *buffer, int n, int omega, int p);


/* In-place FFT Algorithm 1 (Geddes, Czapor, Labahn section 4.7)
 * We require that:
 *  - n is a power of two,
 *  - p is prime,
 *  - coefficients has length n,
 *  - w is a precomputed array of powers of some primitive n-th root of unity
 *    (computed by primitive_root_powers) in Z/p.
 *
 * The result is stored in coefficients; no new memory is allocated.
 */
void fft1(int *const coefficients, int n, const int *const w, int p, int depth);


/* In-place FFT Algorithm 2 (Gerhard and von zur Gathen section 8.2)
 * We require that:
 *  - n is a power of two,
 *  - p is prime,
 *  - coefficients has length n,
 *  - w is a precomputed array of powers of some primitive n-th root of unity
 *    (computed by primitive_root_powers) in Z/p.
 *
 * The result is stored in coefficients; no new memory is allocated.
 */
void fft2(int *const coefficients, int n, const int *const w, int p, int depth);

#ifdef __cplusplus
    }
#endif
