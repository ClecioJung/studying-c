#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"

bool check_if_is_eig(const Matrix A, const double eig, const Vector vec) {
    for (size_t i = 0; i < A.rows; i++) {
        matrix_dec(A, i, i, eig);
    }
    const bool check = matrix_is_null_space(A, vec);
    for (size_t i = 0; i < A.rows; i++) {
        matrix_inc(A, i, i, eig);
    }
    return check;
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
    {
        Vector greatest_vec = (Vector){0};
        const double greatest_eig = matrix_power_method(A, &greatest_vec);
        printf("Greatest eigenvalue: %lg\n", greatest_eig);
        printf("Eigenvector:\n");
        vector_print(greatest_vec);
        if (check_if_is_eig(A, greatest_eig, greatest_vec)) {
            printf("The greatest eigenvalue and eigenvector were properly calculated!\n\n");
        } else {
            fprintf(stderr, "The greatest eigenvalue and eigenvector were NOT properly calculated!\n");
            return EXIT_FAILURE;
        }
        vector_dealloc(&greatest_vec);
    }
    {
        Vector lowest_vec = (Vector){0};
        const double lowest_eig = matrix_inverse_power_method(A, &lowest_vec);
        printf("Lowest eigenvalue: %lg\n", lowest_eig);
        printf("Eigenvector:\n");
        vector_print(lowest_vec);
        if (check_if_is_eig(A, lowest_eig, lowest_vec)) {
            printf("The lowest eigenvalue and eigenvector were properly calculated!\n\n");
        } else {
            fprintf(stderr, "The lowest eigenvalue and eigenvector were NOT properly calculated!\n");
            return EXIT_FAILURE;
        }
        vector_dealloc(&lowest_vec);
    }
    matrix_dealloc(&A);
    return EXIT_SUCCESS;
}
