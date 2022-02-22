#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"
#include "../../lib/scalar.h"

int main(void) {
    Matrix A = matrix_alloc(3, 3);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 4.0);
    matrix_set(A, 0, 1, 12.0);
    matrix_set(A, 0, 2, -16.0);
    matrix_set(A, 1, 0, 12.0);
    matrix_set(A, 1, 1, 37.0);
    matrix_set(A, 1, 2, -43.0);
    matrix_set(A, 2, 0, -16.0);
    matrix_set(A, 2, 1, -43.0);
    matrix_set(A, 2, 2, 98.0);
    matrix_print(A);
    Matrix L = matrix_cholesky_decomposition(A);
    printf("Matrix R:\n");
    matrix_print(L);
    if (matrix_is_lower_triangular(L)) {
        printf("The matrix R is lower-triangular!\n\n");
    } else {
        fprintf(stderr, "The decomposition was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    printf("R * R^T =\n");
    Matrix trR = matrix_transpose(L);
    Matrix mul = matrix_mul(L, trR);
    matrix_print(mul);
    if (matrix_are_equal(mul, A)) {
        printf("This equals to the A matrix!\nSo, we calculated the Cholesky decomposition corectly!\n");
    } else {
        fprintf(stderr, "The decomposition was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    matrix_dealloc(&A);
    matrix_dealloc(&L);
    matrix_dealloc(&trR);
    matrix_dealloc(&mul);
    return EXIT_SUCCESS;
}
