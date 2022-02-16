#include <math.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    Vector v = alloc_vector(2);
    v.data[0] = 1.0;
    v.data[1] = 2.0;
    printf("Vector v:\n");
    print_vector(v);
    Matrix H = householder_matrix(v);
    printf("Matrix H:\n");
    print_matrix(H);
    if (matrix_is_orthogonal(H)) {
        printf("The matrix H is orthogonal!\n");
    }
    free_matrix(&H);
    free_vector(&v);
    return EXIT_SUCCESS;
}
