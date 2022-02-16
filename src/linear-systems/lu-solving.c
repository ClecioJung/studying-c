#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = matrix_alloc(4, 4);
    Vector b = vector_alloc(4);
    printf("Matrix A:\n");
    A.data[0][0] = 10.0;
    A.data[0][1] = -2.0;
    A.data[0][2] = -1.0;
    A.data[0][3] = -1.0;
    A.data[1][0] = -2.0;
    A.data[1][1] = 10.0;
    A.data[1][2] = -1.0;
    A.data[1][3] = -1.0;
    A.data[2][0] = -1.0;
    A.data[2][1] = -1.0;
    A.data[2][2] = 10.0;
    A.data[2][3] = -2.0;
    A.data[3][0] = -1.0;
    A.data[3][1] = -1.0;
    A.data[3][2] = -2.0;
    A.data[3][3] = 10.0;
    matrix_print(A);
    printf("Vector b:\n");
    b.data[0] = 3.0;
    b.data[1] = 15.0;
    b.data[2] = 27.0;
    b.data[3] = -9.0;
    vector_print(b);
    printf("Solution:\n");
    Vector x = lu_solving(A, b);
    vector_print(x);
    vector_dealloc(&x);
    vector_dealloc(&b);
    matrix_dealloc(&A);
    return EXIT_SUCCESS;
}
