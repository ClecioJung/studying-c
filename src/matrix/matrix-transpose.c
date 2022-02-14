#include <math.h>
#include <stdlib.h>

#define MATRIX_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = alloc_matrix(3, 4);
    printf("Matrix A: \n");
    A.data[0][0] = 9.0;
    A.data[0][1] = -18.0;
    A.data[0][2] = 10.0;
    A.data[0][3] = 1.0;
    A.data[1][0] = -36.0;
    A.data[1][1] = 96.0;
    A.data[1][2] = -60.0;
    A.data[1][3] = -5.0;
    A.data[2][0] = 30.0;
    A.data[2][1] = -90.0;
    A.data[2][2] = 60.0;
    A.data[2][3] = 1.0;
    print_matrix(A);
    printf("Transpose matrix: \n");
    Matrix transpose = matrix_transpose(A);
    print_matrix(transpose);
    free_matrix(&A);
    free_matrix(&transpose);
    return EXIT_SUCCESS;
}
