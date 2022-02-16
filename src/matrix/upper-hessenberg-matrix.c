#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = matrix_alloc(4, 4);
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
    matrix_print(A);
    Matrix U, H;
    upper_hessenberg_matrix(A, &U, &H);
    printf("Matrix U:\n");
    matrix_print(U);
    printf("Matrix H:\n");
    matrix_print(H);
    if (matrix_is_orthogonal(U)) {
        printf("The matrix U is orthogonal!\n\n");
    }
    printf("U * H * U^T =\n");
    Matrix transposeU = matrix_transpose(U);
    Matrix mul = matrix_mul_three(U, H, transposeU);
    matrix_print(mul);
    if (matrix_are_equal(mul, A)) {
        printf("This equals to the A matrix!\nSo, we calculated the upper Hessenberg matrix corectly!\n");
    }
    matrix_dealloc(&A);
    matrix_dealloc(&U);
    matrix_dealloc(&H);
    matrix_dealloc(&transposeU);
    matrix_dealloc(&mul);
    return EXIT_SUCCESS;
}
