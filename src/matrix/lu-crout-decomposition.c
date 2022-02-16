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
    Matrix L, U;
    lu_crout_decomposition(A, &L, &U);
    printf("Matrix L:\n");
    print_matrix(L);
    printf("Matrix U:\n");
    print_matrix(U);
    printf("L * U =\n");
    Matrix mul = mul_matrices(L, U);
    print_matrix(mul);
    if (matrices_are_equal(mul, A)) {
        printf("This equals to the A matrix!\nSo, we calculated the LU decomposition corectly!\n");
    }
    free_matrix(&A);
    free_matrix(&L);
    free_matrix(&U);
    free_matrix(&mul);
    return EXIT_SUCCESS;
}
