#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"

int main(void) {
    Vector v = vector_alloc(2);
    v.data[0] = 1.0;
    v.data[1] = 2.0;
    printf("Vector v:\n");
    vector_print(v);
    Matrix H = householder_matrix(v);
    printf("Matrix H:\n");
    matrix_print(H);
    if (matrix_is_orthogonal(H)) {
        printf("The matrix H is orthogonal!\n");
    }
    matrix_dealloc(&H);
    vector_dealloc(&v);
    return EXIT_SUCCESS;
}
