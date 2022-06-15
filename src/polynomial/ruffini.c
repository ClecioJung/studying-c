#include <stdio.h>
#include <stdlib.h>

#include "../../lib/polynomial.h"
#include "../../lib/vector.h"

int main(void) {
    // p(x) = 2*x^3 + 3*x^2 - 4
    Vector p = vector_alloc(4);
    p.data[0] = -4.0;
    p.data[1] = 0.0;
    p.data[2] = 3.0;
    p.data[3] = 2.0;
    // divided by q(x) = x + 1
    Vector q = vector_alloc(2);
    q.data[0] = 1.0;
    q.data[1] = 1.0;
    Vector r = (Vector){0};
    const double remainder = polynomial_ruffini(p, -q.data[0], &r);
    printf("\nUsing Ruffini's rule to divide the polynomial\n");
    polynomial_print(p, "p", 'x');
    printf("  by\n");
    polynomial_print(q, "q", 'x');
    printf("  results in\n");
    polynomial_print(r, "r", 'x');
    printf("  with remainder %lg.\n\n", remainder);
    Vector check = polynomial_multiply(r, q);
    check.data[0] += remainder;
    if (polynomial_are_equal(p, check)) {
        printf("This is the correct result!\n\n");
    } else {
        fprintf(stderr, "The polynomial division was NOT properly calculated\n\n");
        return EXIT_FAILURE;
    }
    vector_dealloc(&p);
    vector_dealloc(&q);
    vector_dealloc(&r);
    vector_dealloc(&check);
    return EXIT_SUCCESS;
}
