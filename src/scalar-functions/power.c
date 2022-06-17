#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define PRECISION 1e-10

typedef struct {
    const double base;
    const uint64_t expoent;
    const double expected;
} Test_Values;

int main(void) {
    const Test_Values values[] = {
        {0.0, 0, 1.0},
        {10.0, 0, 1.0},
        {0.0, 10, 0.0},
        {-1.0, 1, -1.0},
        {-1.0, 2, 1.0},
        {1.0, 1, 1.0},
        {0.5, 2, 0.25},
        {2.0, 2, 4.0},
        {2.0, 3, 8.0},
        {2.0, 4, 16.0},
        {2.0, 5, 32.0},
        {2.0, 10, 1024.0},
        {3.0, 4, 81.0},
        {3.0, 6, 729.0},
        {10.0, 5, 100000.0},
    };
    printf("Power:\n");
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        const double result = power(values[i].base, values[i].expoent);
        if (are_close(result, values[i].expected, PRECISION)) {
            printf("power(%lg, %ld) = %lg\n", values[i].base, values[i].expoent, result);
        } else {
            fprintf(stderr, "Error: imprecise result. power(%lg, %ld) = %lg\n", values[i].base, values[i].expoent, result);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}