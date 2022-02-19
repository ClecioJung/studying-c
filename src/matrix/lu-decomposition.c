#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"
#include "../../lib/scalar.h"

#define PRECISION 1e-10

int main(void) {
    Matrix A = matrix_alloc(5, 5);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 1.0);
    matrix_set(A, 0, 1, 1.0);
    matrix_set(A, 0, 2, 1.0);
    matrix_set(A, 0, 3, 1.0);
    matrix_set(A, 0, 4, 1.0);
    matrix_set(A, 1, 0, 16.0);
    matrix_set(A, 1, 1, 8.0);
    matrix_set(A, 1, 2, 4.0);
    matrix_set(A, 1, 3, 2.0);
    matrix_set(A, 1, 4, 1.0);
    matrix_set(A, 2, 0, 81.0);
    matrix_set(A, 2, 1, 27.0);
    matrix_set(A, 2, 2, 9.0);
    matrix_set(A, 2, 3, 3.0);
    matrix_set(A, 2, 4, 1.0);
    matrix_set(A, 3, 0, 256.0);
    matrix_set(A, 3, 1, 64.0);
    matrix_set(A, 3, 2, 16.0);
    matrix_set(A, 3, 3, 4.0);
    matrix_set(A, 3, 4, 1.0);
    matrix_set(A, 4, 0, 625.0);
    matrix_set(A, 4, 1, 125.0);
    matrix_set(A, 4, 2, 25.0);
    matrix_set(A, 4, 3, 5.0);
    matrix_set(A, 4, 4, 1.0);
    matrix_print(A);
    {
        printf("LU decomposition:\n\n");
        Matrix L, U;
        matrix_lu_decomposition(A, &L, &U);
        printf("Matrix L:\n");
        matrix_print(L);
        if (matrix_is_lower_triangular(L)) {
            printf("The matrix L is lower-triangular!\n\n");
        } else {
            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
            return EXIT_FAILURE;
        }
        printf("Matrix U:\n");
        matrix_print(U);
        if (matrix_is_upper_triangular(U)) {
            printf("The matrix U is upper-triangular!\n\n");
        } else {
            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
            return EXIT_FAILURE;
        }
        printf("L * U =\n");
        Matrix mul = matrix_mul(L, U);
        matrix_print(mul);
        if (matrix_are_equal(mul, A)) {
            printf("This equals to the A matrix!\nSo, we calculated the LU decomposition corectly!\n\n");
        } else {
            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
            return EXIT_FAILURE;
        }
        {
            printf("LU decomposition (version with less memory usage):\n");
            printf("In this version, the LU matrices are stored in the matrix A, changing its contents!\n\n");
            Matrix copy = matrix_copy(A);
            matrix_lu_dec_over(copy);
            printf("Matrix L/U:\n");
            matrix_print(copy);
            for (size_t i = 0; i < copy.rows; i++) {
                for (size_t j = 0; j < copy.cols; j++) {
                    if (j >= i) {
                        if (!are_close(matrix_get(copy, i, j), matrix_get(U, i, j), PRECISION)) {
                            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
                            return EXIT_FAILURE;
                        }
                    } else {
                        if (!are_close(matrix_get(copy, i, j), matrix_get(L, i, j), PRECISION)) {
                            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
            matrix_undo_lu_over(copy);
            printf("Matrix A recovered:\n");
            matrix_print(copy);
            if (!matrix_are_equal(copy, A)) {
                fprintf(stderr, "The decomposition was NOT properly calculated!\n");
                return EXIT_FAILURE;
            }
            matrix_dealloc(&copy);
        }
        matrix_dealloc(&L);
        matrix_dealloc(&U);
        matrix_dealloc(&mul);
    }
    {
        printf("LU Crout decomposition:\n\n");
        Matrix L, U;
        matrix_lu_crout_decomposition(A, &L, &U);
        printf("Matrix L:\n");
        matrix_print(L);
        if (matrix_is_lower_triangular(L)) {
            printf("The matrix L is lower-triangular!\n\n");
        } else {
            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
            return EXIT_FAILURE;
        }
        printf("Matrix U:\n");
        matrix_print(U);
        if (matrix_is_upper_triangular(U)) {
            printf("The matrix U is upper-triangular!\n\n");
        } else {
            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
            return EXIT_FAILURE;
        }
        printf("L * U =\n");
        Matrix mul = matrix_mul(L, U);
        matrix_print(mul);
        if (matrix_are_equal(mul, A)) {
            printf("This equals to the A matrix!\nSo, we calculated the LU Crout decomposition corectly!\n\n");
        } else {
            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
            return EXIT_FAILURE;
        }
        {
            printf("LU Crout decomposition (version with less memory usage):\n");
            printf("In this version, the LU matrices are stored in the matrix A, changing its contents!\n\n");
            Matrix copy = matrix_copy(A);
            matrix_lu_crout_dec_over(copy);
            printf("Matrix L/U:\n");
            matrix_print(copy);
            for (size_t i = 0; i < copy.rows; i++) {
                for (size_t j = 0; j < copy.cols; j++) {
                    if (j > i) {
                        if (!are_close(matrix_get(copy, i, j), matrix_get(U, i, j), PRECISION)) {
                            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
                            return EXIT_FAILURE;
                        }
                    } else {
                        if (!are_close(matrix_get(copy, i, j), matrix_get(L, i, j), PRECISION)) {
                            fprintf(stderr, "The decomposition was NOT properly calculated!\n");
                            return EXIT_FAILURE;
                        }
                    }
                }
            }
            matrix_undo_lu_crout_over(copy);
            printf("Matrix A recovered:\n");
            matrix_print(copy);
            if (!matrix_are_equal(copy, A)) {
                fprintf(stderr, "The decomposition was NOT properly calculated!\n");
                return EXIT_FAILURE;
            }
            matrix_dealloc(&copy);
        }
        matrix_dealloc(&L);
        matrix_dealloc(&U);
        matrix_dealloc(&mul);
    }
    matrix_dealloc(&A);
    return EXIT_SUCCESS;
}
