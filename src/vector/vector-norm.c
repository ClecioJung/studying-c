#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"
#include "../../lib/vector.h"

#define PRECISION 1e-10

int main(void) {
    const double expected_result = 5.0;
    printf("Vector a:\n");
    Vector a = vector_new(3, 0.0, 3.0, 4.0);
    vector_print(a);
    const double result = vector_norm(a);
    if (are_close(result, expected_result, PRECISION)) {
        printf("Euclidean norm: %lg\n", result);
    } else {
        fprintf(stderr, "The euclidean norm was NOT properly calculated: %lg\n", result);
        return EXIT_FAILURE;
    }
    vector_dealloc(&a);
    return EXIT_SUCCESS;
}
