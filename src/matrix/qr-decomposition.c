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
    Matrix Q, R;
    qr_decomposition(A, &Q, &R);
    printf("Matrix Q:\n");
    matrix_print(Q);
    printf("Matrix R:\n");
    matrix_print(R);
    if (matrix_is_orthogonal(Q)) {
        printf("The matrix Q is orthogonal!\n\n");
    }
    printf("Q * R =\n");
    Matrix mul = matrix_mul(Q, R);
    matrix_print(mul);
    if (matrix_are_equal(mul, A)) {
        printf("This equals to the A matrix!\nSo, we calculated the QR decomposition corectly!\n");
    }
    matrix_dealloc(&A);
    matrix_dealloc(&Q);
    matrix_dealloc(&R);
    matrix_dealloc(&mul);
    return EXIT_SUCCESS;
}
