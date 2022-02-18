#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define EULER 2.7182818284590452354
#define MAX_ITERATIONS 10000
#define PRECISION 1e-10

double euler(void) {
    double e = 1.0;
    for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
        const double delta = 1.0 / ((double)factorial(i));
        e += delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return e;
}

int main(void) {
    const double result = euler();
    if (are_close(result, EULER, PRECISION)) {
        printf("Euler constant was aproximated to be %lg\n", result);
        return EXIT_SUCCESS;
    } else {
        fprintf(stderr, "Couldn't approximate euler's constant: %lg\n", result);
        return EXIT_FAILURE;
    }
}
