#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = alloc_matrix(2, 3);
    printf("Matrix A:\n");
    A.data[0][0] = 1.0;
    A.data[0][1] = 2.0;
    A.data[0][2] = 3.0;
    A.data[1][0] = 4.0;
    A.data[1][1] = 5.0;
    A.data[1][2] = 6.0;
    print_matrix(A);
    printf("Pseudo inverse matrix:\n");
    Matrix pseudo_inv = pseudo_inverse(A);
    print_matrix(pseudo_inv);
    printf("pseudo_inverse(A) * A =\n");
    Matrix mul_left = mul_matrices(pseudo_inv, A);
    print_matrix(mul_left);
    printf("A * pseudo_inverse(A) =\n");
    Matrix mul_right = mul_matrices(A, pseudo_inv);
    print_matrix(mul_right);
    printf("A * pseudo_inverse(A) * A =\n");
    Matrix mul = mul_3_matrices(A, pseudo_inv, A);
    print_matrix(mul);
    if (matrices_are_equal(mul, A)) {
        printf("This equals to the A matrix!\nSo, we calculated the pseudo-inverse matrix corectly!\n");
    }
    free_matrix(&mul);
    free_matrix(&A);
    free_matrix(&pseudo_inv);
    free_matrix(&mul_left);
    free_matrix(&mul_right);
    return EXIT_SUCCESS;
}
