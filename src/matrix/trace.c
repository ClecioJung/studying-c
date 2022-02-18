#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"
#include "../../lib/scalar.h"

#define PRECISION 1e-10

int main(void) {
    const double expected_result = 7.0;
    Matrix A = matrix_alloc(3, 3);
    printf("Matrix A:\n");
    matrix_set(A, 0, 0, 2.0);
    matrix_set(A, 0, 1, 3.0);
    matrix_set(A, 0, 2, -1.0);
    matrix_set(A, 1, 0, 4.0);
    matrix_set(A, 1, 1, 4.0);
    matrix_set(A, 1, 2, -3.0);
    matrix_set(A, 2, 0, 2.0);
    matrix_set(A, 2, 1, -3.0);
    matrix_set(A, 2, 2, 1.0);
    matrix_print(A);
    const double result = trace(A);
    if (are_close(result, expected_result, PRECISION)) {
        printf("Trace: %lg\n", result);
    } else {
        fprintf(stderr, "The trace was NOT properly calculated: %lg\n", result);
        return EXIT_FAILURE;
    }
    matrix_dealloc(&A);
    return EXIT_SUCCESS;
}
