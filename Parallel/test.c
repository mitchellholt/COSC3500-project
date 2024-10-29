#include "fin_field.h"
#include <stdio.h>
#include <emmintrin.h>


int main(void) {
    __m128i vec1 = _mm_set_epi32(40694, 13612, 19386, 2347);
    __m128i vec2 = _mm_set_epi32(26723, 42533, 27999, 5860);
    __m128i vec3 = _mm_set_epi32(61273, 56058, 38161, 50355);

    __m128i res;
    int *arr = (int*)&res;

    // Should be [32120, 61199, 59647, 64095]
    //__m128i res = vec_triple_mulmod(vec1, vec2, vec3, 65537);

    for (int i = 0; i < 4; i++) {
        printf("%d ", arr[4 - i - 1]);
    }
    printf("\n");

    // try the other computation as well
    res = vec_mulmod(vec_mulmod(vec1, vec2, 65537), vec3, 65537);
    for (int i = 0; i < 4; i++) {
        printf("%d ", arr[4 - i - 1]);
    }
    printf("\n");

    return 0;
}
