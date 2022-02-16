#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
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
    printf("Dot product: %g\n", dot_product(a, b));
    vector_dealloc(&a);
    vector_dealloc(&b);
    return EXIT_SUCCESS;
}
