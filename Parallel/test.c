#include "fin_field.h"
#include <stdio.h>
#include <emmintrin.h>


int main(void) {
    int coefficients[] = { 1, 2, 3, 4 };
    __m128i coeffs = _mm_loadu_si128((__m128i*)coefficients);
    int ww[] = { 1, 8, 64, 512 };
    const int p = 65537;

    // Correct method
    int s1 = coefficients[0];
    int s2 = coefficients[1];
    int t1 = mulmod(ww[0], coefficients[2], p);
    int t2 = mulmod(ww[1], coefficients[3], p);
    coefficients[0] = addmod(s1, t1, p);
    coefficients[1] = addmod(s2, t2, p);
    coefficients[2] = submod(s1, t1, p);
    coefficients[3] = submod(s2, t2, p);

    for (int i = 0; i < 4; i++) {
        printf(" %d", coefficients[i]);
    }
    printf("\n");

    // Scuffed method
    __m128i w = _mm_shuffle_epi32(_mm_loadu_si128((__m128i*)ww), _MM_SHUFFLE(1, 0, 1, 0));
    __m128i t = _mm_shuffle_epi32(coeffs, _MM_SHUFFLE(3, 2, 3, 2));
    __m128i s = _mm_shuffle_epi32(coeffs, _MM_SHUFFLE(1, 0, 1, 0));
    t = vec_mulmod(t, w, p);
    __m128i adds = vec_addmod(s, t, p);
    __m128i subs = vec_submod(s, t, p);
    coeffs = _mm_blend_epi32(adds, subs, 0xC); // 0b1100

    int *arr = (int*)&coeffs;
    for (int i = 0; i < 4; i++) {
        printf(" %d", arr[i]);
    }
    printf("\n");

    return 0;
}
