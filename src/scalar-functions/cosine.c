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
        {-3.14159265358979323846, -1.0},
        {-2.356194490192345, -0.7071067811865475},
        {-1.5707963267948966, 0.0},
        {-0.7853981633974483, 0.7071067811865475},
        {0.0, 1.0},
        {0.7853981633974483, 0.7071067811865475},
        {1.5707963267948966, 0.0},
        {2.356194490192345, -0.7071067811865475},
        {3.14159265358979323846, -1.0},
        {100*3.14159265358979323846, 1.0},
    };
    printf("Cosine:\n");
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        const double result = cosine(values[i].input);
        if (are_close(result, values[i].expected, PRECISION)) {
            printf("cosine(%lg) = %lg\n", values[i].input, result);
        } else {
            fprintf(stderr, "Error: imprecise result. cosine(%lg) = %lg\n", values[i].input, result);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}