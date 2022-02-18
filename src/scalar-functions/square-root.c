#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define PRECISION 1e-10

typedef struct {
    const double input;
    const double expected;
} Test_Values;

int main(void) {
    const Test_Values values[] = {
        {0.0, 0},
        {1.0, 1},
        {2.0, 1.4142135623730951},
        {3.0, 1.7320508075688772},
        {4.0, 2},
        {5.0, 2.23606797749979},
        {6.0, 2.449489742783178},
        {7.0, 2.6457513110645907},
        {8.0, 2.8284271247461903},
        {9.0, 3},
        {10.0, 3.1622776601683795},
    };
    printf("Square root:\n");
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        const double result = square_root(values[i].input);
        if (are_close(result, values[i].expected, PRECISION)) {
            printf("square_root(%lg) = %lg\n", values[i].input, result);
        } else {
            fprintf(stderr, "Error: imprecise result. square_root(%lg) = %lg\n", values[i].input, result);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}