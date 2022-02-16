#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Vector a = alloc_vector(3);
    Vector b = alloc_vector(3);
    printf("Vector a:\n");
    a.data[0] = 2.0;
    a.data[1] = 3.0;
    a.data[2] = 4.0;
    print_vector(a);
    printf("Vector b:\n");
    b.data[0] = 5.0;
    b.data[1] = 6.0;
    b.data[2] = 7.0;
    print_vector(b);
    printf("Cross product:\n");
    Vector cross = cross_product(a, b);
    print_vector(cross);
    free_vector(&a);
    free_vector(&b);
    free_vector(&cross);
    return EXIT_SUCCESS;
}
