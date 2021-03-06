#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"

int main(void) {
    Matrix A = matrix_alloc(3, 3);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 9.0);
    matrix_set(A, 0, 1, -18.0);
    matrix_set(A, 0, 2, 10.0);
    matrix_set(A, 1, 0, -36.0);
    matrix_set(A, 1, 1, 96.0);
    matrix_set(A, 1, 2, -60.0);
    matrix_set(A, 2, 0, 30.0);
    matrix_set(A, 2, 1, -90.0);
    matrix_set(A, 2, 2, 60.0);
    matrix_print(A);
    printf("Inverse matrix:\n");
    Matrix inv = matrix_inverse(A);
    matrix_print(inv);
    printf("A * matrix_inverse(A) =\n");
    Matrix mul = matrix_mul(A, inv);
    matrix_print(mul);
    Matrix identity = matrix_identity(A.rows);
    if (matrix_are_equal(mul, identity)) {
        printf("The inverse matrix was properly calculated!\n");
    } else {
        fprintf(stderr, "The inverse matrix was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    matrix_dealloc(&identity);
    matrix_dealloc(&A);
    matrix_dealloc(&inv);
    matrix_dealloc(&mul);
    return EXIT_SUCCESS;
}
