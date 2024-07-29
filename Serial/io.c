#include <stdio.h>
#include "io.h"


void print_int_array(FILE *stream, const int *const arr, int len) {
    fprintf(stream, "[");
    for (int i = 0; i < len - 1; i++) {
        fprintf(stream, "%d, ", arr[i]);
    }
    fprintf(stream, "%d]", arr[len - 1]);
}
