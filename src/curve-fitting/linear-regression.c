#include <stdio.h>
#include <stdlib.h>

#include "../../lib/curve-fitting.h"
#include "../../lib/scalar.h"

#define PRECISION 1e-10

int main(void) {
    Vector x = vector_alloc(3);
    Vector y = vector_alloc(3);
    x.data[0] = 0.0;
    x.data[1] = 1.0;
    x.data[2] = 2.0;
    y.data[0] = 5.0;
    y.data[1] = 7.0;
    y.data[2] = 9.0;
    double a, b;
    printf("Linear regression resulted in:\n");
    double r2 = linear_regression(x, y, &a, &b);
    printf("a: %lg\n", a);
    printf("b: %lg\n", b);
    printf("The coefficient of determination is (r^2): %lg\n", r2);
    for (size_t i = 0; i < x.len; i++) {
        const double result = a * x.data[i] + b;
        if (are_close(result, y.data[i], PRECISION)) {
            printf("The approx. polynomial avaliated at %lg results in %lg, while the original value was %lg\n", x.data[i], result, y.data[i]);
        } else {
            fprintf(stderr, "Linear regression: couldn't interpolate correctly the provided values!\n");
            return EXIT_FAILURE;
        }
    }
    vector_dealloc(&x);
    vector_dealloc(&y);
    return EXIT_SUCCESS;
}