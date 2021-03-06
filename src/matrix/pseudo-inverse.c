#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"

int main(void) {
    Matrix A = matrix_alloc(2, 3);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 1.0);
    matrix_set(A, 0, 1, 2.0);
    matrix_set(A, 0, 2, 3.0);
    matrix_set(A, 1, 0, 4.0);
    matrix_set(A, 1, 1, 5.0);
    matrix_set(A, 1, 2, 6.0);
    matrix_print(A);
    printf("Pseudo inverse matrix:\n");
    Matrix pseudo_inv = matrix_pseudo_inverse(A);
    matrix_print(pseudo_inv);
    printf("pseudo_inverse(A) * A =\n");
    Matrix mul_left = matrix_mul(pseudo_inv, A);
    matrix_print(mul_left);
    printf("A * pseudo_inverse(A) =\n");
    Matrix mul_right = matrix_mul(A, pseudo_inv);
    matrix_print(mul_right);
    printf("A * pseudo_inverse(A) * A =\n");
    Matrix mul = matrix_mul_three(A, pseudo_inv, A);
    matrix_print(mul);
    if (matrix_are_equal(mul, A)) {
        printf("This equals to the A matrix!\nSo, he pseudo-inverse matrix was properly calculated!\n");
    } else {
        fprintf(stderr, "The pseudo-inverse matrix was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    matrix_dealloc(&mul);
    matrix_dealloc(&A);
    matrix_dealloc(&pseudo_inv);
    matrix_dealloc(&mul_left);
    matrix_dealloc(&mul_right);
    return EXIT_SUCCESS;
}
