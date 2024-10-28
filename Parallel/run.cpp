#include <stdint.h>
#include <mm_malloc.h>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <chrono>
#include "polynomial.h"
#include "io.h"
#include "fin_field.h"

#define BAD_ARGS_ERR 1
#define BAD_FILE_ERR 2
#define FRICK        3 // something real bad happened

#define USAGE_STRING "Usage: run <num_multiplications_per_size> [-o filename] {powers}\n"
#define NUM_MULTS_INDEX    1
#define OUT_INDEX          2
#define MIN_ARGS           3

#define ALIGN        64
#define FERMAT_PRIME 65537
#define MAX_POWER_2  16

// Example usage: ./run 10 out 10 11 (runs 10 examples of sizes 2^10 and 2^10)
int main(int argc, char **argv) {
    int num_mults;
    int powers_start_index;
    FILE *out = NULL;

    if (argc < MIN_ARGS || sscanf(argv[NUM_MULTS_INDEX], "%d", &num_mults) != 1) {
        fprintf(stderr, USAGE_STRING);
        return BAD_ARGS_ERR;
    }

    if (strcmp(argv[OUT_INDEX], "-o") == 0) {
        out = fopen(argv[OUT_INDEX + 1], "w");
        if (!out) {
            perror("Opening provided file for writing");
            return BAD_FILE_ERR;
        }
        powers_start_index = OUT_INDEX + 2;
    } else {
        powers_start_index = OUT_INDEX;
    }

    // Read input powers and check they are correct sizes
    int num_sizes = argc - powers_start_index;
    int *sizes = (int*)_mm_malloc(sizeof(int) * num_sizes, ALIGN);
    int biggest = 0;
    for (int i = 0; i < num_sizes; i++) {
        if (sscanf(argv[powers_start_index + i], "%d", sizes + i) != 1
                || sizes[i] < 8
                || sizes[i] > 16) {
            fprintf(stderr, "argv[%d] is not between 8 and 16\n",
                    powers_start_index + i);
            _mm_free(sizes);
            if (out) fclose(out);
            return BAD_ARGS_ERR;
        }
        sizes[i] = (int)(1 << (uint32_t)sizes[i]);
        if (sizes[i] > biggest) biggest = sizes[i];
    }

    // Allocate buffers
    int *a_buf = (int*)_mm_malloc(biggest * sizeof(int), ALIGN);
    int *b_buf = (int*)_mm_malloc(biggest * sizeof(int), ALIGN);

    srand(time(NULL));

    for (int i = 0; i < num_sizes; i++) {
        int n = sizes[i];
        int n2 = n/2;
        int omega = fermat_primitive_root(n, FERMAT_PRIME);
        if (!omega) {
            if (out) fclose(out);
            _mm_free(a_buf);
            _mm_free(b_buf);
            _mm_free(sizes);

            fprintf(
                stderr,
                "BUG couldn't create a primitive root of unity when n is %d\n",
                n);
            return FRICK;
        }

        // record the total number of microseconds spent in 
        long totalTime = 0;

        for (int j = 0; j < num_mults; j++) {
            // create random polys and print to file
            rand_poly(a_buf, n2, FERMAT_PRIME);
            rand_poly(b_buf, n2, FERMAT_PRIME);
            memset(a_buf + n2, 0, n2 * sizeof(int));
            memset(b_buf + n2, 0, n2 * sizeof(int));

            if (out) { // DEBUG mode
                print_int_array(out, a_buf, n2);
                fprintf(out, "\n");
                print_int_array(out, b_buf, n2);
                fprintf(out, "\n");
            }

            int expected_degree = degree(n, a_buf) + degree(n, b_buf);

            // multiply them and print result
            const auto startTime = std::chrono::high_resolution_clock::now();
            fast_multiply(a_buf, b_buf, omega, n, FERMAT_PRIME);
            const auto endTime = std::chrono::high_resolution_clock::now();
            totalTime += std::chrono::duration_cast<std::chrono::microseconds>(endTime - startTime).count();

            if (out) { // DEBUG mode
                print_int_array(out, a_buf, n);
                fprintf(out, "\n");
            }

            int recovered_degree = degree(n, a_buf);
            if (recovered_degree != expected_degree) {
                fprintf(stderr,
                        "Error: expected a product of degree %d, but got %d\n",
                        expected_degree,
                        recovered_degree);
            }
        }

        printf("%d multiplications of size %d. Total time %ld us\n",
                num_mults, n, totalTime);
    }

    if (out) fclose(out);
    _mm_free(a_buf);
    _mm_free(b_buf);
    _mm_free(sizes);
}
