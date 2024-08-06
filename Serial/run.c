#include <stdint.h>
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
#define FRICK        3

#define FERMAT_PRIME 65537
#define MAX_POWER_2  16


// Nice struct to make it easier to reuse memory
struct buffer {
    int capacity; // number of ints allocated for data
    int *data;       // actual data
};


static inline void buffer_realloc(struct buffer *buf, int n) {
    if (buf->capacity >= n) return;
    buf->data = realloc(buf->data, sizeof(int) * n);
    buf->capacity = n;
}


static inline uint32_t pow2(uint32_t k) {
    return 1 << k;
}


int main(int argc, char **argv) {
    int num_mults;

    if (argc != 3 || sscanf(argv[1], "%d", &num_mults) != 1) {
        fprintf(stderr, "Usage: run <num_multiplications> <output file_name>\n");
        return BAD_ARGS_ERR;
    }

    FILE *out = fopen(argv[2], "w");
    if (!out) {
        perror("Opening provided file for writing");
        return BAD_FILE_ERR;
    }

    struct buffer a_buf = { 0, NULL };
    struct buffer b_buf = { 0, NULL };

    srand(time(NULL));

    for (int i = 0; i < num_mults; i++) {
        //int n = (int)pow2( ((uint32_t)rand() % 8) + 8 ); // between 2^9 and 2^16
        int n = (int)pow2(16);

        // allocate buffers (if necessary)
        buffer_realloc(&a_buf, n);
        buffer_realloc(&b_buf, n);

        // create random polys and print to file
        int n2 = n/2;
        rand_poly(a_buf.data, n2, FERMAT_PRIME);
        rand_poly(b_buf.data, n2, FERMAT_PRIME);
        memset(a_buf.data + n2, 0, n2);
        memset(b_buf.data + n2, 0, n2);

        print_int_array(out, a_buf.data, n2);
        fprintf(out, "\n");
        print_int_array(out, b_buf.data, n2);
        fprintf(out, "\n");

        // multiply them and print result
        int omega = fermat_primitive_root(n, FERMAT_PRIME);
        if (!omega) {
            fclose(out);
            free(a_buf.data);
            free(b_buf.data);

            fprintf(
                stderr,
                "FRICK couldn't create a primitive root of unity when n is %d\n",
                n);
            return FRICK;
        }

        fast_multiply(a_buf.data, b_buf.data, omega, n, FERMAT_PRIME);
        print_int_array(out, a_buf.data, n);
        fprintf(out, "\n");
    }

    fclose(out);
    free(a_buf.data);
    free(b_buf.data);
}
