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
    Matrix P, D;
    matrix_diagonalization(A, &P, &D);
    printf("Matrix D:\n");
    matrix_print(D);
    if (matrix_is_diagonal(D)) {
        printf("The matrix D is diagonal!\n\n");
    } else {
        fprintf(stderr, "The diagonalization was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    printf("Matrix P:\n");
    matrix_print(P);
    printf("P * D * inv(P) =\n");
    Matrix Pinv = matrix_inverse(P);
    Matrix mul = matrix_mul_three(P, D, Pinv);
    matrix_print(mul);
    if (matrix_are_equal(mul, A)) {
        printf("This equals to the A matrix!\nSo, we calculated the diagonalization corectly!\n");
    } else {
        fprintf(stderr, "The diagonalization was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    matrix_dealloc(&A);
    matrix_dealloc(&P);
    matrix_dealloc(&D);
    matrix_dealloc(&Pinv);
    matrix_dealloc(&mul);
    return EXIT_SUCCESS;
}
