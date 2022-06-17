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
        {-10.0, -0.6483608274590866},
        {-5.0, 3.380515006246586},
        {-3.14159265358979323846, 0.0},
        {-0.7853981633974483, -1.0},
        {0.0, 0.0},
        {0.7853981633974483, 1.0},
        {3.14159265358979323846, 0.0},
        {100*3.14159265358979323846, 0.0},
    };
    printf("Tangent:\n");
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        const double result = tangent(values[i].input);
        if (are_close(result, values[i].expected, PRECISION)) {
            printf("tangent(%lg) = %lg\n", values[i].input, result);
        } else {
            fprintf(stderr, "Error: imprecise result. tangent(%lg) = %lg\n", values[i].input, result);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}