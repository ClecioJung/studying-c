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
    {
        Vector p = vector_alloc(4);
        p.data[0] = -2.0;
        p.data[1] = -1.0;
        p.data[2] = 0.0;
        p.data[3] = 2.0;
        if (test_algorithm(p) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        vector_dealloc(&p);
    }
    {
        Vector p = vector_alloc(4);
        p.data[0] = -1.0;
        p.data[1] = 3.0;
        p.data[2] = -3.0;
        p.data[3] = 1.0;
        if (test_algorithm(p) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        vector_dealloc(&p);
    }
    {
        Vector p = vector_alloc(5);
        p.data[0] = 6.0;
        p.data[1] = -17.0;
        p.data[2] = 17.0;
        p.data[3] = -7.0;
        p.data[4] = 1.0;
        if (test_algorithm(p) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        vector_dealloc(&p);
    }
    {
        Vector p = vector_alloc(5);
        p.data[0] = -8.0;
        p.data[1] = 4.0;
        p.data[2] = 6.0;
        p.data[3] = -5.0;
        p.data[4] = 1.0;
        if (test_algorithm(p) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        vector_dealloc(&p);
    }
    {
        Vector p = vector_alloc(5);
        p.data[0] = 0.001;
        p.data[1] = -1.0111;
        p.data[2] = 11.1111;
        p.data[3] = -11.101;
        p.data[4] = 1.0;
        if (test_algorithm(p) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        vector_dealloc(&p);
    }
    {
        Vector p = vector_alloc(9);
        p.data[0] = 0.531441;
        p.data[1] = -2.480058;
        p.data[2] = 3.287061;
        p.data[3] = 1.56006;
        p.data[4] = -7.1685;
        p.data[5] = 4.32;
        p.data[6] = 2.35;
        p.data[7] = -3.4;
        p.data[8] = 1.0;
        if (test_algorithm(p) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        vector_dealloc(&p);
    }
    {
        Vector p = vector_alloc(11);
        p.data[0] = 0.0;
        p.data[1] = 0.0;
        p.data[2] = 0.0;
        p.data[3] = -0.9504;
        p.data[4] = 6.7512;
        p.data[5] = -20.5012;
        p.data[6] = 34.5004;
        p.data[7] = -34.75;
        p.data[8] = 20.95;
        p.data[9] = -7.0;
        p.data[10] = 1.0;
        if (test_algorithm(p) == EXIT_FAILURE) {
            return EXIT_FAILURE;
        }
        vector_dealloc(&p);
    }
    return EXIT_SUCCESS;
}
