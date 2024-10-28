#include "fin_field.h"
#include <stdio.h>
#include <emmintrin.h>


int main(void) {
    __m128i vec1 = _mm_set_epi32(13612, 19386, 2347, 26723);
    __m128i vec2 = _mm_set_epi32(42533, 27999, 5860, 61273);

    // Should be [5338, 11180, 56187, 21971]
    __m128i res = vec_mulmod(vec1, vec2, 65537);
    int *arr = (int*)&res;

    for (int i = 0; i < 4; i++) {
        printf("%d ", arr[4 - i - 1]);
    }
    printf("\n");

    return 0;
}
