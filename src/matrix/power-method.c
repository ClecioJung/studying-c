#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"

int main(void) {
    Matrix A = matrix_alloc(4, 4);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 1.0);
    matrix_set(A, 0, 1, 1.0);
    matrix_set(A, 0, 2, 1.0);
    matrix_set(A, 0, 3, 1.0);
    matrix_set(A, 1, 0, 8.0);
    matrix_set(A, 1, 1, 4.0);
    matrix_set(A, 1, 2, 2.0);
    matrix_set(A, 1, 3, 1.0);
    matrix_set(A, 2, 0, 27.0);
    matrix_set(A, 2, 1, 9.0);
    matrix_set(A, 2, 2, 3.0);
    matrix_set(A, 2, 3, 1.0);
    matrix_set(A, 3, 0, 64.0);
    matrix_set(A, 3, 1, 16.0);
    matrix_set(A, 3, 2, 4.0);
    matrix_set(A, 3, 3, 1.0);
    matrix_print(A);
    Vector vec = (Vector){0};
    const double eig = power_method(A, &vec);
    printf("Greatest eigenvalue: %lg\n", eig);
    printf("Eigenvector:\n");
    vector_print(vec);
    {
        for (size_t i = 0; i < A.rows; i++) {
            matrix_dec(A, i, i, eig);
        }
        if (matrix_is_null_space(A, vec)) {
            printf("The eigenvalue and eigenvector were properly calculated!\n\n");
        } else {
            fprintf(stderr, "The eigenvalue and eigenvector were NOT properly calculated!\n");
            return EXIT_FAILURE;
        }
    }
    matrix_dealloc(&A);
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}
