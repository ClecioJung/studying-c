#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"

Vector eigenvector(const Matrix A, const double eigenvalue) {
    if (!matrix_is_squared(A)) {
        return (Vector){0};  // Invalid operation
    }
    Matrix A_copy = matrix_copy(A);
    Vector zeros = vector_init(A.rows, 0.0);
    for (size_t i = 0; i < A.rows; i++) {
        matrix_dec(A_copy, i, i, eigenvalue);
    }
    Matrix pinv = matrix_pseudo_inverse(A_copy);
    Vector eigenvector = matrix_mul_vector(pinv, zeros);
    matrix_dealloc(&A_copy);
    matrix_dealloc(&pinv);
    vector_dealloc(&zeros);
    return eigenvector;
}

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
    printf("Eigenvalues:\n");
    Vector eig = matrix_eigenvalues(A);
    vector_print(eig);
    for (size_t i = 0; i < eig.len; i++) {
        Vector vec = eigenvector(A, eig.data[i]);
        vector_print(vec);
        vector_dealloc(&vec);
    }
    matrix_dealloc(&A);
    vector_dealloc(&eig);
    return EXIT_SUCCESS;
}
