#include <stdio.h>
#include <stdlib.h>

#include "../../lib/curve-fitting.h"

int main(void) {
    Vector x = vector_alloc(3);
    Vector y = vector_alloc(3);
    x.data[0] = 0.0;
    x.data[1] = 1.0;
    x.data[2] = 2.0;
    y.data[0] = 0.0;
    y.data[1] = 2.0;
    y.data[2] = 4.0;
    const double value = 1.5;
    printf("Lagrange interpolation at %g results in %g\n", value, lagrange_interpolation(x, y, value));
    vector_dealloc(&x);
    vector_dealloc(&y);
    return EXIT_SUCCESS;
}
