#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Vector x = alloc_vector(5);
    Vector y = alloc_vector(5);
    x.data[0] = 0.0;
    x.data[1] = 1.0;
    x.data[2] = 2.0;
    x.data[3] = 3.0;
    x.data[4] = 4.0;
    y.data[0] = 5.05;
    y.data[1] = 6.9;
    y.data[2] = 9.1;
    y.data[3] = 12.3;
    y.data[4] = 15.1;
    printf("Polynomial regression resulted in: \n");
    Vector pol = polynomial_regression(x, y, 2);
    print_vector(pol);
    for (size_t i = 0; i < x.len; i++) {
        printf("The approx. polynomial avaliated at %g results in %g, while the original vlaue was %g\n", x.data[i], compute_polynomial(pol, x.data[i]), y.data[i]);
    }

    free_vector(&pol);
    free_vector(&x);
    free_vector(&y);
    return EXIT_SUCCESS;
}