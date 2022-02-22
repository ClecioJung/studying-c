#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "../../lib/linear-systems.h"
#include "../../lib/scalar.h"

#define PRECISION 1e-9

typedef System_Type (*Linear_Systems_Fn)(const Matrix A, const Vector b, Vector *const x);

typedef struct {
    const char *const name;
    const Linear_Systems_Fn method_fn;
} Linear_Systems_Algorithm;

int test_algorithm(const Matrix A, const Vector b, const Linear_Systems_Algorithm algorithm) {
    struct timeval stop, start;
    Vector x;
    gettimeofday(&start, NULL);
    algorithm.method_fn(A, b, &x);
    gettimeofday(&stop, NULL);
    Vector result = matrix_mul_vector(A, x);
    bool correctness = vector_are_equal(result, b);
    if (correctness) {
        unsigned long int delta_us = (stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec);
        printf("Took %04ld us for the %s to find the correct solution!\n", delta_us, algorithm.name);
    } else {
        fprintf(stderr, "The %s didn't found the correct solution!\n", algorithm.name);
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
    };
    const size_t len = 100;
    const double min = 0.0;
    const double max = 100.0;
    Matrix A = matrix_random(len, len, min, max);
    Vector b = vector_random(len, min, max);
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(algorithms[0]); i++) {
        if (test_algorithm(A, b, algorithms[i]) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    vector_dealloc(&b);
    matrix_dealloc(&A);
    return EXIT_SUCCESS;
}
