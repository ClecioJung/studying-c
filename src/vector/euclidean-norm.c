#include <math.h>
#include <stdlib.h>

#define MATRIX_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Vector a = alloc_vector(3);
    printf("Vector a: \n");
    a.data[0] = 4.0;
    a.data[1] = 3.0;
    a.data[2] = 0.0;
    print_vector(a);
    printf("Euclidean norm: %g\n", euclidean_norm(a));
    free_vector(&a);
    return EXIT_SUCCESS;
}
