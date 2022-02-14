#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = alloc_matrix(3, 3);
    printf("Matrix A: \n");
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
    printf("Inverse matrix: \n");
    Matrix inv = matrix_inverse(A);
    print_matrix(inv);
    printf("A * matrix_inverse(A) =\n");
    Matrix mul = mul_matrices(A, inv);
    print_matrix(mul);
    Matrix identity = identity_matrix(A.rows);
    if (matrices_are_equal(mul, identity)) {
        printf("This equals to the identity matrix!\nSo, we calculated the inverse matrix corectly!\n");
    }
    free_matrix(&identity);
    free_matrix(&A);
    free_matrix(&inv);
    free_matrix(&mul);
    return EXIT_SUCCESS;
}
