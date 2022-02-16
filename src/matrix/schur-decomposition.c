#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = alloc_matrix(4, 4);
    printf("Matrix A:\n");
    A.data[0][0] = 1.0;
    A.data[0][1] = 1.0;
    A.data[0][2] = 1.0;
    A.data[0][3] = 1.0;
    A.data[1][0] = 8.0;
    A.data[1][1] = 4.0;
    A.data[1][2] = 2.0;
    A.data[1][3] = 1.0;
    A.data[2][0] = 27.0;
    A.data[2][1] = 9.0;
    A.data[2][2] = 3.0;
    A.data[2][3] = 1.0;
    A.data[3][0] = 64.0;
    A.data[3][1] = 16.0;
    A.data[3][2] = 4.0;
    A.data[3][3] = 1.0;
    print_matrix(A);
    Matrix U, T;
    schur_decomposition(A, &U, &T);
    printf("Matrix U:\n");
    print_matrix(U);
    printf("Matrix T:\n");
    print_matrix(T);
    if (matrix_is_orthogonal(U)) {
        printf("The matrix U is orthogonal!\n\n");
    }
    printf("U * T * U^T =\n");
    Matrix transposeU = matrix_transpose(U);
    Matrix mul = mul_3_matrices(U, T, transposeU);
    print_matrix(mul);
    if (matrices_are_equal(mul, A)) {
        printf("This equals to the A matrix!\nSo, we calculated the Schur decomposition corectly!\n");
    }
    free_matrix(&A);
    free_matrix(&U);
    free_matrix(&T);
    free_matrix(&transposeU);
    free_matrix(&mul);
    return EXIT_SUCCESS;
}
