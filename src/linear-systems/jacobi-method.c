#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = alloc_matrix(4, 4);
    Vector b = alloc_vector(4);
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
    print_matrix(A);
    printf("Vector b:\n");
    b.data[0] = 3.0;
    b.data[1] = 15.0;
    b.data[2] = 27.0;
    b.data[3] = -9.0;
    print_vector(b);
    printf("Columns condition: %d\n", columns_condition(A));
    printf("Rows condition: %d\n", rows_condition(A));
    printf("Sassenfeld condition: %d\n\n", sassenfeld_condition(A));
    printf("Resoution by Jacobi method:\n");
    Vector x = jacobi_method(A, b);
    print_vector(x);
    free_vector(&x);
    free_vector(&b);
    free_matrix(&A);
    return EXIT_SUCCESS;
}
