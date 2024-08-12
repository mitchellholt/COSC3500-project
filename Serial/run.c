#include <stdint.h>
#include <mm_malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "polynomial.h"
#include "io.h"
#include "fin_field.h"

#define BAD_ARGS_ERR 1
#define BAD_FILE_ERR 2
#define FRICK        3 // something real bad happened

#define USAGE_STRING "Usage: run <num_multiplications_per_size> <output file_name> {powers}\n"
#define NUM_MULTS_INDEX    1
#define OUT_FILE_INDEX     2
#define POWERS_START_INDEX 3
#define MIN_ARGS           4

#define ALIGN        64
#define FERMAT_PRIME 65537
#define MAX_POWER_2  16


static inline uint32_t pow2(uint32_t k) {
    return 1 << k;
}


int main(int argc, char **argv) {
    int num_mults;

    if (argc < MIN_ARGS || sscanf(argv[NUM_MULTS_INDEX], "%d", &num_mults) != 1) {
        fprintf(stderr, USAGE_STRING);
        return BAD_ARGS_ERR;
    }

    FILE *out = fopen(argv[OUT_FILE_INDEX], "w");
    if (!out) {
        perror("Opening provided file for writing");
        return BAD_FILE_ERR;
    }

    // Read input powers and check they are correct sizes
    int num_sizes = argc - POWERS_START_INDEX;
    int *sizes = _mm_malloc(sizeof(int) * num_sizes, ALIGN);
    int biggest = 0;
    for (int i = 0; i < num_sizes; i++) {
        if (sscanf(argv[POWERS_START_INDEX + i], "%d", sizes + i) != 1
                || sizes[i] < 1
                || sizes[i] > 16) {
            fprintf(stderr, "argv[%d] is not between 1 and 16\n",
                    POWERS_START_INDEX + i);
            free(sizes);
            fclose(out);
            return BAD_ARGS_ERR;
        }
        sizes[i] = (int)pow2((uint32_t)sizes[i]);
        if (sizes[i] > biggest) biggest = sizes[i];
    }

    // Allocate buffers
    int *a_buf = _mm_malloc(biggest * sizeof(int), ALIGN);
    int *b_buf = _mm_malloc(biggest * sizeof(int), ALIGN);

    srand(time(NULL));

    for (int i = 0; i < num_sizes; i++) {
        int n = sizes[i];
        int n2 = n/2;
        int omega = fermat_primitive_root(n, FERMAT_PRIME);
        if (!omega) {
            fclose(out);
            free(a_buf);
            free(b_buf);
            free(sizes);

            fprintf(
                stderr,
                "BUG couldn't create a primitive root of unity when n is %d\n",
                n);
            return FRICK;
        }
        fprintf(out, "%d\n%d\n", n, num_mults); // used for verification

        for (int j = 0; j < num_mults; j++) {
            // create random polys and print to file
            rand_poly(a_buf, n2, FERMAT_PRIME);
            rand_poly(b_buf, n2, FERMAT_PRIME);
            memset(a_buf + n2, 0, n2 * sizeof(int));
            memset(b_buf + n2, 0, n2 * sizeof(int));

            print_int_array(out, a_buf, n2);
            fprintf(out, "\n");
            print_int_array(out, b_buf, n2);
            fprintf(out, "\n");

            // multiply them and print result
            fast_multiply(a_buf, b_buf, omega, n, FERMAT_PRIME);
            print_int_array(out, a_buf, n);
            fprintf(out, "\n");
        }
    }

    fclose(out);
    free(a_buf);
    free(b_buf);
    free(sizes);
}
