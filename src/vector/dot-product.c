#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"
#include "../../lib/vector.h"

#define PRECISION 1e-10

int main(void) {
    const double expected_result = 4.0;
    Vector a = vector_alloc(3);
    Vector b = vector_alloc(3);
    printf("Vector a:\n");
    a.data[0] = 1.0;
    a.data[1] = 2.0;
    a.data[2] = 3.0;
    vector_print(a);
    printf("Vector b:\n");
    b.data[0] = 2.0;
    b.data[1] = 1.0;
    b.data[2] = 0.0;
    vector_print(b);
    const double result = dot_product(a, b);
    if (are_close(result, expected_result, PRECISION)) {
        printf("Dot product: %lg\n", result);
    } else {
        fprintf(stderr, "The dot product was NOT properly calculated: %lg\n", result);
        return EXIT_FAILURE;
    }
    vector_dealloc(&a);
    vector_dealloc(&b);
    return EXIT_SUCCESS;
}
