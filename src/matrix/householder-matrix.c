#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/matrix.h"

int main(void) {
    Vector v = vector_alloc(5);
    v.data[0] = 1.0;
    v.data[1] = 2.0;
    v.data[2] = 3.0;
    v.data[3] = 4.0;
    v.data[4] = 5.0;
    printf("Vector v:\n");
    vector_print(v);
    Matrix H = householder_matrix(v);
    printf("Matrix H:\n");
    matrix_print(H);
    if (matrix_is_orthogonal(H) && matrix_is_symmetric(H)) {
        printf("The Householder matrix is orthogonal and symmetric!\n");
    } else {
        fprintf(stderr, "The Householder matrix was NOT properly calculated\n");
        return EXIT_FAILURE;
    }
    matrix_dealloc(&H);
    vector_dealloc(&v);
    return EXIT_SUCCESS;
}
