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
    printf("Symmetric matrix:\n");
    Matrix sym = matrix_symmetric(A);
    matrix_print(sym);
    printf("Skew-symmetric matrix:\n");
    Matrix skew = matrix_skew_symmetric(A);
    matrix_print(skew);
    printf("Sum of the previous matrices:\n");
    Matrix sum = matrix_sum(sym, skew);
    matrix_print(sum);
    if (matrix_is_symmetric(sym) && matrix_is_skew_symmetric(skew) && matrix_are_equal(sum, A)) {
        printf("The symmetric and skew-symmetric matrices were properly calculated!\n");
    } else {
        fprintf(stderr, "The symmetric and skew-symmetric matrices were NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    matrix_dealloc(&A);
    matrix_dealloc(&sym);
    matrix_dealloc(&skew);
    matrix_dealloc(&sum);
    return EXIT_SUCCESS;
}
