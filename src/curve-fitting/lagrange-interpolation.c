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
    y.data[0] = 0.0;
    y.data[1] = 2.0;
    y.data[2] = 4.0;
    const data_type value = 1.5;
    printf("Lagrange interpolation at %g results in %g\n", value, lagrange_interpolation(x, y, value));
    free_vector(&x);
    free_vector(&y);
    return EXIT_SUCCESS;
}
