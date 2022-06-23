#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"
#include "../../lib/scalar.h"

// Remember to free the returned vector after calling this function!
Vector f(const Vector x) {
    if (x.len != 3) {
        return (Vector){0};  // Invalid operation
    }
    Vector y = vector_alloc(4);
    if (vector_is_valid(y)) {
        y.data[0] = x.data[0];
        y.data[1] = 5.0 * x.data[2];
        y.data[2] = 4.0 * x.data[1] * x.data[1] - 2.0 * x.data[2];
        y.data[3] = x.data[2] * sine(x.data[0]);
    }
    return y;
}

int main(void) {
    Vector x = vector_alloc(3);
    x.data[0] = 0.0;
    x.data[1] = 1.0;
    x.data[2] = -1.0;
    Matrix J = matrix_jacobian(f, x);
    printf("The Jacobian matrix is:\n");
    matrix_print(J);
    vector_dealloc(&x);
    matrix_dealloc(&J);
    return EXIT_SUCCESS;
}
