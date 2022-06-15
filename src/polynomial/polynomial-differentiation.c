#include <stdio.h>
#include <stdlib.h>

#include "../../lib/polynomial.h"
#include "../../lib/vector.h"

int main(void) {
    // p(x) = 2x4 + 3x – 2
    Vector p = vector_alloc(5);
    p.data[0] = -2.0;
    p.data[1] = 3.0;
    p.data[2] = 0.0;
    p.data[3] = 0.0;
    p.data[4] = 2.0;
    const double x = 2.0;
    printf("\nConsidering the polynomial\n");
    polynomial_print(p, "p", 'x');
    printf("\nIt's value at %lg is %lg\n", x, polynomial_evaluation(p, x));
    printf("It's first derivative evaluated at %lg is %lg\n", x, polynomial_first_diff(p, x));
    printf("It's second derivative evaluated at %lg is %lg\n", x, polynomial_diff(p, 2, x));
    printf("It's third derivative evaluated at %lg is %lg\n", x, polynomial_diff(p, 3, x));
    printf("It's fourth derivative evaluated at %lg is %lg\n", x, polynomial_diff(p, 4, x));
    printf("It's fifth derivative evaluated at %lg is %lg\n\n", x, polynomial_diff(p, 5, x));
    vector_dealloc(&p);
    return EXIT_SUCCESS;
}
