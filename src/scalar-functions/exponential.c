#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

int main(void) {
    printf("Exponential:\n");
    for (int i = -9; i < 10; i++) {
        printf("exponential(%02d) = %g\n", i, exponential((double)i));
    }
    return EXIT_SUCCESS;
}
