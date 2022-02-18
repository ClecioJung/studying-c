#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"
#include "../../lib/vector.h"

#define PRECISION 1e-10

int main(void) {
    Vector a = vector_alloc(3);
    Vector b = vector_alloc(3);
    printf("Vector a:\n");
    a.data[0] = 2.0;
    a.data[1] = 3.0;
    a.data[2] = 4.0;
    vector_print(a);
    printf("Vector b:\n");
    b.data[0] = 5.0;
    b.data[1] = 6.0;
    b.data[2] = 7.0;
    vector_print(b);
    printf("Cross product:\n");
    Vector cross = cross_product(a, b);
    vector_print(cross);
    if (vector_are_orthogonal(a, cross) && vector_are_orthogonal(a, cross)) {
        printf("The cross product is orthogonal to the vectors a and b!\n");
    } else {
        fprintf(stderr, "The cross product was NOT properly calculated!\n");
        return EXIT_FAILURE;
    }
    vector_dealloc(&a);
    vector_dealloc(&b);
    vector_dealloc(&cross);
    return EXIT_SUCCESS;
}
