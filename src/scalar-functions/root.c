#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define PRECISION 1e-10

typedef struct {
    const double input;
    const uint16_t n;
    const double expected;
} Test_Values;

int main(void) {
    const Test_Values values[] = {
        {100.0, 0, 1.0},
        {0.0, 1, 0.0},
        {1.0, 2, 1.0},
        {2.0, 2, 1.4142135623730951},
        {3.0, 2, 1.7320508075688772},
        {2.0, 3, 1.2599210498948732},
        {3.0, 3, 1.4422495703074083},
        {8.0, 3, 2.0},
        {81.0, 4, 3.0},
        {97.65625, 5, 2.5},
        {1024.0, 10, 2.0},
    };
    printf("Root:\n");
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        const double result = root(values[i].input, values[i].n);
        if (are_close(result, values[i].expected, PRECISION)) {
            printf("root(%lg, %d) = %lg\n", values[i].input, values[i].n, result);
        } else {
            fprintf(stderr, "Error: imprecise result. root(%lg, %d) = %lg\n", values[i].input, values[i].n, result);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}