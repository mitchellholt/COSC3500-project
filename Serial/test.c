#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "fin_field.h"
#include "fft.h"
#include "io.h"
#include "polynomial.h"

#define BAD_ARGS_ERR 1


int main(int argc, char **argv) {
    if (argc < 2) {
        fprintf(stderr, "usage: test [primitive|fft|polymult] ...\n");
        return BAD_ARGS_ERR;
    }

    int n,p,omega;

    if (!strcmp(argv[1], "primitive")) {
        if (argc != 4) {
            fprintf(stderr, "Usage: %s primitive <n> <p>\n", argv[0]);
            return BAD_ARGS_ERR;
        }
        n = atoi(argv[1]);
        p = atoi(argv[2]);

        omega = primitive_root(n, p);
        printf("omega: %d\n", omega);

        int *w = primitive_root_powers(n, omega, p);
        print_int_array(stdout, w, n);
        printf("\n");

    } else if (!strcmp(argv[1], "fft")) {
        n = 16;
        p = 97;
        omega = primitive_root(n, p);
        printf("omega: %d\n", omega);
        int *w = primitive_root_powers(n, omega, p);
        print_int_array(stdout, w, n);
        printf("\n");
        {
            // a(x) from MATH895 assignment 1
            int coefficients[16] = { 1, 39, 0, 57, 11, 0, 19, 0, 1, 0, 0, 0, 0, 0, 0, 0 };
            fft1(coefficients, n, w, p);
            printf("fft1: ");
            print_int_array(stdout, coefficients, 16);
            printf("\n");
        }
        {
            int coefficients[16] = { 1, 39, 0, 57, 11, 0, 19, 0, 1, 0, 0, 0, 0, 0, 0, 0 };
            fft2(coefficients, n, w, p);
            printf("fft2: ");
            print_int_array(stdout, coefficients, 16);
            printf("\n");
        }

    } else if (!strcmp(argv[1], "polymult")) {
        // Example: Q1 from MATH895 A1
        {
            int a[16] = { 1, 39, 0, 57, 11, 0, 19, 0, 1, 0, 0, 0, 0, 0, 0, 0 };
            int b[16] = { 7, 0,  0, 22, 44, 0, 17, 9, 0, 0, 0, 0, 0, 0, 0, 0 };
            if (!fast_multiply(a, b, 16, 97)) {
                fprintf(stderr, "Couldn't find primitive root of unity\n");
                return 2;
            }
            print_int_array(stdout, a, 16);
            printf("\n");
        }

        // Example 4.8 from Geddes, Czapor, Labahn
        {
            int a[8] = { 1,  37, 1, 3, 0, 0, 0, 0 };
            int b[8] = { 38, 5,  2, 1, 0, 0, 0, 0 };
            if (!fast_multiply(a, b, 8, 41)) {
                fprintf(stderr, "Couldn't find primitive root of unity\n");
                return 2;
            }
            print_int_array(stdout, a, 8);
            printf("\n");
        }
    } else {
        fprintf(stderr, "usage: test [primitive|fft|polymult] ...\n");
        return BAD_ARGS_ERR;
    }

    return 0;
}
