#include <stdio.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Vector x = alloc_vector(3);
    Vector y = alloc_vector(3);
    x.data[0] = 0.0;
    x.data[1] = 1.0;
    x.data[2] = 2.0;
    y.data[0] = 5.05;
    y.data[1] = 6.9;
    y.data[2] = 9.1;
    data_type a, b;
    printf("Linear regression resulted in:\n");
    data_type r = linear_regression(x, y, &a, &b);
    printf("a: %g\n", a);
    printf("b: %g\n", b);
    printf("correlation coefficient (r): %g\n", r);
    free_vector(&x);
    free_vector(&y);
    return EXIT_SUCCESS;
}