#include "../../lib/linear-systems.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define PRECISION 1e-9

typedef Vector (*Linear_Systems_Fn)(const Matrix A, const Vector b);

typedef struct {
    const char *const name;
    const Linear_Systems_Fn method_fn;
} Linear_Systems_Algorithm;

int test_algorithm(const Matrix A, const Vector b, const Linear_Systems_Algorithm algorithm) {
    Vector x = algorithm.method_fn(A, b);
    printf("The vector solution found by the %s is:\n", algorithm.name);
    vector_print(x);
    Vector result = matrix_mul_vector(A, x);
    bool correctness = vector_are_equal(result, b);
    if (correctness) {
        printf("This is the correct solution!\n");
    } else {
        fprintf(stderr, "This is not the correct solution!\n");
    }
    vector_dealloc(&x);
    vector_dealloc(&result);
    return correctness ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(void) {
    const Linear_Systems_Algorithm algorithms[] = {
        {"Gaussian elimination", gaussian_elimination},
        {"Gauss-Jordan method", gauss_jordan},
        {"LU method", lu_solving},
        {"Jacobi method", jacobi_method},
        {"Gauss-Seidel method", gauss_seidel},
    };
    Matrix A = matrix_alloc(4, 4);
    Vector b = vector_alloc(4);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 10.0);
    matrix_set(A, 0, 1, -2.0);
    matrix_set(A, 0, 2, -1.0);
    matrix_set(A, 0, 3, -1.0);
    matrix_set(A, 1, 0, -2.0);
    matrix_set(A, 1, 1, 10.0);
    matrix_set(A, 1, 2, -1.0);
    matrix_set(A, 1, 3, -1.0);
    matrix_set(A, 2, 0, -1.0);
    matrix_set(A, 2, 1, -1.0);
    matrix_set(A, 2, 2, 10.0);
    matrix_set(A, 2, 3, -2.0);
    matrix_set(A, 3, 0, -1.0);
    matrix_set(A, 3, 1, -1.0);
    matrix_set(A, 3, 2, -2.0);
    matrix_set(A, 3, 3, 10.0);
    matrix_print(A);
    printf("Vector b:\n");
    b.data[0] = 3.0;
    b.data[1] = 15.0;
    b.data[2] = 27.0;
    b.data[3] = -9.0;
    vector_print(b);
    {
        bool condition = columns_condition(A);
        if (condition) {
            printf("Columns condition returned true!\n");
        } else {
            fprintf(stderr, "Expected the condition to be true!");
            return EXIT_FAILURE;
        }
        condition = rows_condition(A);
        if (condition) {
            printf("Rows condition returned true!\n");
        } else {
            fprintf(stderr, "Expected the condition to be true!");
            return EXIT_FAILURE;
        }
        condition = sassenfeld_condition(A);
        if (condition) {
            printf("Sassenfeld condition returned true!\n\n");
        } else {
            fprintf(stderr, "Expected the condition to be true!");
            return EXIT_FAILURE;
        }
    }
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(algorithms[0]); i++) {
        if (test_algorithm(A, b, algorithms[i]) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    vector_dealloc(&b);
    matrix_dealloc(&A);
    return EXIT_SUCCESS;
}
