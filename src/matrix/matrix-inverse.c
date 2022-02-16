#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = matrix_alloc(3, 3);
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
    matrix_print(A);
    printf("Inverse matrix:\n");
    Matrix inv = matrix_inverse(A);
    matrix_print(inv);
    printf("A * matrix_inverse(A) =\n");
    Matrix mul = matrix_mul(A, inv);
    matrix_print(mul);
    Matrix identity = matrix_identity(A.rows);
    if (matrix_are_equal(mul, identity)) {
        printf("This equals to the identity matrix!\nSo, we calculated the inverse matrix corectly!\n");
    }
    matrix_dealloc(&identity);
    matrix_dealloc(&A);
    matrix_dealloc(&inv);
    matrix_dealloc(&mul);
    return EXIT_SUCCESS;
}
