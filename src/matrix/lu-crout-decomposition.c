#include <math.h>
#include <stdlib.h>

#define MATRIX_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = alloc_matrix(4, 4);
    printf("Matrix A: \n");
    A.data[0][0] = 2.0;
    A.data[0][1] = 4.0;
    A.data[0][2] = 1.0;
    A.data[0][3] = -3.0;
    A.data[1][0] = -1.0;
    A.data[1][1] = -1.0;
    A.data[1][2] = 2.0;
    A.data[1][3] = 4.0;
    A.data[2][0] = 4.0;
    A.data[2][1] = 2.0;
    A.data[2][2] = -3.0;
    A.data[2][3] = 5.0;
    A.data[3][0] = 5.0;
    A.data[3][1] = -4.0;
    A.data[3][2] = -3.0;
    A.data[3][3] = 1.0;
    print_matrix(A);
    Matrix L, U;
    lu_crout_decomposition(A, &L, &U);
    printf("Matrix L: \n");
    print_matrix(L);
    printf("Matrix U: \n");
    print_matrix(U);
    printf("L * U =\n");
    Matrix mul = mul_matrices(L, U);
    print_matrix(mul);
    free_matrix(&A);
    free_matrix(&L);
    free_matrix(&U);
    free_matrix(&mul);
    return EXIT_SUCCESS;
}
