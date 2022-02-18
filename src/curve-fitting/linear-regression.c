#include <stdio.h>
#include <stdlib.h>

#include "../../lib/curve-fitting.h"

int main(void) {
    Vector x = vector_alloc(3);
    Vector y = vector_alloc(3);
    x.data[0] = 0.0;
    x.data[1] = 1.0;
    x.data[2] = 2.0;
    y.data[0] = 5.05;
    y.data[1] = 6.9;
    y.data[2] = 9.1;
    double a, b;
    printf("Linear regression resulted in:\n");
    double r = linear_regression(x, y, &a, &b);
    printf("a: %g\n", a);
    printf("b: %g\n", b);
    printf("correlation coefficient (r): %g\n", r);
    vector_dealloc(&x);
    vector_dealloc(&y);
    return EXIT_SUCCESS;
}