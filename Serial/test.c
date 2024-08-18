#include <stdio.h>
#include <stdbool.h>
#include "fin_field.h"
#include "fft.h"
#include "io.h"
#include "polynomial.h"


int main(int argc, char **argv) {
    if (argc == 1) {
        // Example: Q1 from MATH895 A1
        {
            int a[16] = { 1, 39, 0, 57, 11, 0, 19, 0, 1, 0, 0, 0, 0, 0, 0, 0 };
            int b[16] = { 7, 0,  0, 22, 44, 0, 17, 9, 0, 0, 0, 0, 0, 0, 0, 0 };
            int omega = 8;

            fast_multiply(a, b, omega, 16, 97);
            print_int_array(stdout, a, 16);
            printf("\n");
        }

        // Example 4.8 from Geddes, Czapor, Labahn
        {
            int a[8] = { 1,  37, 1, 3, 0, 0, 0, 0 };
            int b[8] = { 38, 5,  2, 1, 0, 0, 0, 0 };
            int omega = 14;

            fast_multiply(a, b, omega, 8, 41);
            print_int_array(stdout, a, 8);
            printf("\n");
        }
    } else {
        printf("Not implemented\n");
        return 0;
    }

    return 0;
}
