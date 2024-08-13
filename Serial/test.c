#include <stdio.h>
#include <stdint.h>
#include "fin_field.h"


int main(int argc, char **argv) {
    int a,b;
    if (argc != 3 || sscanf(argv[1], "%d", &a) != 1 || sscanf(argv[1], "%d", &b) != 1) {
        fprintf(stderr, "Usage: test <num> <num>\n");
        return 1;
    }

    int64_t t = (a * b) % FERMAT_PRIME;
    int expected = (int) t;
    int computed = mulmod(a, b);
    
    if (expected != computed) {
        printf("Expected %d, but got %d\n", expected, computed);
    } else {
        printf("success");
    }
    return 0;
}
