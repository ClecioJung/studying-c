#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"

int main(void) {
    Matrix A = matrix_alloc(3, 4);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 9.0);
    matrix_set(A, 0, 1, -18.0);
    matrix_set(A, 0, 2, 10.0);
    matrix_set(A, 0, 3, 1.0);
    matrix_set(A, 1, 0, -36.0);
    matrix_set(A, 1, 1, 96.0);
    matrix_set(A, 1, 2, -60.0);
    matrix_set(A, 1, 3, -5.0);
    matrix_set(A, 2, 0, 30.0);
    matrix_set(A, 2, 1, -90.0);
    matrix_set(A, 2, 2, 60.0);
    matrix_set(A, 2, 3, 1.0);
    matrix_print(A);
    printf("Transpose matrix:\n");
    Matrix transpose = matrix_transpose(A);
    matrix_print(transpose);
    matrix_dealloc(&A);
    matrix_dealloc(&transpose);
    return EXIT_SUCCESS;
}
