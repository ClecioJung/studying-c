#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Vector a = alloc_vector(3);
    Vector b = alloc_vector(3);
    printf("Vector a:\n");
    a.data[0] = 1.0;
    a.data[1] = 2.0;
    a.data[2] = 3.0;
    print_vector(a);
    printf("Vector b:\n");
    b.data[0] = 2.0;
    b.data[1] = 1.0;
    b.data[2] = 0.0;
    print_vector(b);
    printf("Dot product: %g\n", dot_product(a, b));
    free_vector(&a);
    free_vector(&b);
    return EXIT_SUCCESS;
}
