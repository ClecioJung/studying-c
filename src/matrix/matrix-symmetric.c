#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = alloc_matrix(3, 3);
    printf("Matrix A:\n");
    A.data[0][0] = 9.0;
    A.data[0][1] = -18.0;
    A.data[0][2] = 10.0;
    A.data[1][0] = -36.0;
    A.data[1][1] = 96.0;
    A.data[1][2] = -60.0;
    A.data[2][0] = 30.0;
    A.data[2][1] = -90.0;
    A.data[2][2] = 60.0;
    print_matrix(A);
    printf("Symmetric matrix:\n");
    Matrix sym = matrix_symmetric(A);
    print_matrix(sym);
    printf("Skew-symmetric matrix:\n");
    Matrix skew = matrix_skew_symmetric(A);
    print_matrix(skew);
    printf("Sum of the previous matrices:\n");
    Matrix sum = sum_matrices(sym, skew);
    print_matrix(sum);
    if (matrices_are_equal(sum, A)) {
        printf("This equals to the A matrix!\nSo, we calculated the Symmetric and Skew-symmetric matrices corectly!\n");
    }
    free_matrix(&A);
    free_matrix(&sym);
    free_matrix(&skew);
    free_matrix(&sum);
    return EXIT_SUCCESS;
}
