#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MATRIX_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    printf("Exponential:\n");
    for (int i = -9; i < 10; i++) {
        printf("exponential(%02d) = %g\n", i, exponential((data_type)i));
    }
    return EXIT_SUCCESS;
}
