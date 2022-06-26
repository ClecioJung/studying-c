#include <stdio.h>
#include <stdlib.h>

#include "../../lib/polynomial.h"
#include "../../lib/vector.h"

int main(void) {
    for (size_t i = 0; i < 11; i++) {
        printf("The Chebyshev polynomial of the first kind and order %ld is\n", i);
        Vector p = polynomial_chebyshev_first_kind(i);
        polynomial_print(p, "T", 'x');
        vector_dealloc(&p);
    }
    printf("\n");
    for (size_t i = 0; i < 11; i++) {
        printf("The Chebyshev polynomial of the second kind and order %ld is\n", i);
        Vector p = polynomial_chebyshev_second_kind(i);
        polynomial_print(p, "U", 'x');
        vector_dealloc(&p);
    }
    return EXIT_SUCCESS;
}
