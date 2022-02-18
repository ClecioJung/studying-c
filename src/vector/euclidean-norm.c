#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/vector.h"

int main(void) {
    Vector a = vector_alloc(3);
    printf("Vector a:\n");
    a.data[0] = 4.0;
    a.data[1] = 3.0;
    a.data[2] = 0.0;
    vector_print(a);
    printf("Euclidean norm: %g\n", euclidean_norm(a));
    vector_dealloc(&a);
    return EXIT_SUCCESS;
}
