#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/complex-vector.h"
#include "../../lib/polynomial.h"
#include "../../lib/vector.h"

int test_algorithm(const Vector p) {
    printf("\nEvaluating the polynomial\n");
    polynomial_print(p, "p", 'x');
    Complex_Vector roots = polynomial_find_roots(p);
    printf("\nThe roots found are:\n");
    complex_vector_print(roots);
    for (size_t i = 0; i < roots.len; i++) {
        if (!complex_is_null(polynomial_complex_evaluation(p, roots.data[i]))) {
            fprintf(stderr, "The polynomial root was NOT properly calculated\n\n");
            return EXIT_FAILURE;
        }
    }
    printf("This is the correct result!\n");
    complex_vector_dealloc(&roots);
    return EXIT_SUCCESS;
}

int main(void) {
    Vector p = vector_alloc(8);
    p.data[0] = -48.0;
    p.data[1] = 220.0;
    p.data[2] = -416.0;
    p.data[3] = 419.0;
    p.data[4] = -242.0;
    p.data[5] = 80.0;
    p.data[6] = -14.0;
    p.data[7] = 1.0;
    printf("\nEvaluating the polynomial\n");
    polynomial_print(p, "p", 'x');
    Complex_Vector roots = polynomial_find_roots(p);
    for (size_t i = 0; i < roots.len; i++) {
        const uint16_t multiplicity = polynomial_root_multiplicity(p, roots.data[i]);
        printf("The multiplicity of the root ");
        complex_print(roots.data[i]);
        printf(" is %d!\n", multiplicity);
    }
    complex_vector_dealloc(&roots);
    vector_dealloc(&p);
    return EXIT_SUCCESS;
}
