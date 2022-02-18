#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/curve-fitting.h"
#include "../../lib/scalar.h"

#define PRECISION 1e-10

int main(void) {
    Vector x = vector_alloc(5);
    Vector y = vector_alloc(5);
    x.data[0] = 0.0;
    x.data[1] = 1.0;
    x.data[2] = 2.0;
    x.data[3] = 3.0;
    x.data[4] = 4.0;
    y.data[0] = 2.0;
    y.data[1] = 4.0;
    y.data[2] = 8.0;
    y.data[3] = 14.0;
    y.data[4] = 22.0;
    printf("Polynomial regression resulted in:\n");
    Vector pol = polynomial_regression(x, y, 2);
    vector_print(pol);
    for (size_t i = 0; i < x.len; i++) {
        const double result = compute_polynomial(pol, x.data[i]);
        if (are_close(result, y.data[i], PRECISION)) {
            printf("The approx. polynomial avaliated at %lg results in %lg, while the original value was %lg\n", x.data[i], result, y.data[i]);
        } else {
            fprintf(stderr, "Polynomial regression: couldn't interpolate correctly the provided values!\n");
            return EXIT_FAILURE;
        }
    }
    vector_dealloc(&pol);
    vector_dealloc(&x);
    vector_dealloc(&y);
    return EXIT_SUCCESS;
}