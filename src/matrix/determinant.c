#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Matrix A = alloc_matrix(3, 3);
    printf("Matrix A: \n");
    A.data[0][0] = 2.0;
    A.data[0][1] = 3.0;
    A.data[0][2] = -1.0;
    A.data[1][0] = 4.0;
    A.data[1][1] = 4.0;
    A.data[1][2] = -3.0;
    A.data[2][0] = 2.0;
    A.data[2][1] = -3.0;
    A.data[2][2] = 1.0;
    print_matrix(A);
    printf("Determinant: %g\n", determinant(A));
    free_matrix(&A);
    return EXIT_SUCCESS;
}
