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
        {-10.0, -1.4711276743037347},
        {-2.0, -1.1071487177940904},
        {-1.0, -0.7853981633974483},
        {-0.5, -0.4636476090008061},
        {0.0, 0.0},
        {0.5, 0.4636476090008061},
        {1.0, 0.7853981633974483},
        {2.0, 1.1071487177940904},
        {10.0, 1.4711276743037347},
    };
    printf("Arctangent:\n");
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        const double result = arctangent(values[i].input);
        if (are_close(result, values[i].expected, PRECISION)) {
            printf("arctangent(%lg) = %lg\n", values[i].input, result);
        } else {
            fprintf(stderr, "Error: imprecise result. arctangent(%lg) = %lg\n", values[i].input, result);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}