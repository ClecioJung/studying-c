#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"

int main(void) {
    Matrix A = matrix_alloc(4, 4);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 1.0);
    matrix_set(A, 0, 1, 1.0);
    matrix_set(A, 0, 2, 1.0);
    matrix_set(A, 0, 3, 1.0);
    matrix_set(A, 1, 0, 8.0);
    matrix_set(A, 1, 1, 4.0);
    matrix_set(A, 1, 2, 2.0);
    matrix_set(A, 1, 3, 1.0);
    matrix_set(A, 2, 0, 27.0);
    matrix_set(A, 2, 1, 9.0);
    matrix_set(A, 2, 2, 3.0);
    matrix_set(A, 2, 3, 1.0);
    matrix_set(A, 3, 0, 64.0);
    matrix_set(A, 3, 1, 16.0);
    matrix_set(A, 3, 2, 4.0);
    matrix_set(A, 3, 3, 1.0);
    matrix_print(A);
    Matrix U, T;
    schur_decomposition(A, &U, &T);
    printf("Matrix U:\n");
    matrix_print(U);
    if (matrix_is_orthogonal(U)) {
        printf("The matrix U is orthogonal!\n\n");
    } else {
        fprintf(stderr, "The decomposition was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    printf("Matrix T:\n");
    matrix_print(T);
    if (matrix_is_upper_triangular(T)) {
        printf("The matrix T is upper-triangular!\n\n");
    } else {
        fprintf(stderr, "The decomposition was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    printf("U * T * U^T =\n");
    Matrix transposeU = matrix_transpose(U);
    Matrix mul = matrix_mul_three(U, T, transposeU);
    matrix_print(mul);
    if (matrix_are_equal(mul, A)) {
        printf("This equals to the A matrix!\nSo, we calculated the Schur decomposition corectly!\n");
    } else {
        fprintf(stderr, "The decomposition was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    matrix_dealloc(&A);
    matrix_dealloc(&U);
    matrix_dealloc(&T);
    matrix_dealloc(&transposeU);
    matrix_dealloc(&mul);
    return EXIT_SUCCESS;
}
