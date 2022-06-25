#include <stdio.h>
#include <stdlib.h>

#include "../../lib/polynomial.h"
#include "../../lib/vector.h"

int main(void) {
    for (size_t i = 0; i < 11; i++) {
        printf("The Legendre polynomial of order %ld is\n", i);
        Vector p = polynomial_legendre(i);
        polynomial_print(p, "p", 'x');
        printf("with roots:\n");
        Vector roots = polynomial_legendre_roots(i);
        vector_print(roots);
        vector_dealloc(&p);
        vector_dealloc(&roots);
    }
    return EXIT_SUCCESS;
}
