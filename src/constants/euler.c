#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

data_type euler(void) {
    data_type e = 1.0;
    for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
        const data_type delta = 1.0 / ((data_type)factorial(i));
        e += delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return e;
}

int main(void) {
    printf("Euler constant was aproximated to be %lg\n", euler());
    return EXIT_SUCCESS;
}
