#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

int main(void) {
    printf("Square root:\n");
    for (int i = 0; i <= 10; i++) {
        printf("square_root(%d) = %g\n", i, square_root((double)i));
    }
    return EXIT_SUCCESS;
}
