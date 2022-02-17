#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = matrix_alloc(3, 3);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 2.0);
    matrix_set(A, 0, 1, 3.0);
    matrix_set(A, 0, 2, -1.0);
    matrix_set(A, 1, 0, 4.0);
    matrix_set(A, 1, 1, 4.0);
    matrix_set(A, 1, 2, -3.0);
    matrix_set(A, 2, 0, 2.0);
    matrix_set(A, 2, 1, -3.0);
    matrix_set(A, 2, 2, 1.0);
    matrix_print(A);
    printf("Trace: %g\n", trace(A));
    matrix_dealloc(&A);
    return EXIT_SUCCESS;
}
