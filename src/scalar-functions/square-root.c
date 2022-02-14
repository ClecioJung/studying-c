#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MATRIX_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    printf("Square root:\n");
    for (int i = 0; i <= 10; i++) {
        printf("square_root(%d) = %g\n", i, square_root((data_type)i));
    }
    return EXIT_SUCCESS;
}
