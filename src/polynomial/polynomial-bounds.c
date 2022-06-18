#include <stdio.h>
#include <stdlib.h>

#include "../../lib/polynomial.h"
#include "../../lib/vector.h"

int main(void) {
    // p(x) = 3*x^6 + 4*x^3 – 2*x^2 – 6
    Vector p = vector_alloc(7);
    p.data[0] = -6.0;
    p.data[1] = 0.0;
    p.data[2] = -2.0;
    p.data[3] = 4.0;
    p.data[4] = 0.0;
    p.data[5] = 0.0;
    p.data[6] = 3.0;
    printf("\nEvaluating the polynomial\n");
    polynomial_print(p, "p", 'x');
    printf("\nThe Cauchy upper bound is %lg\n", polynomial_cauchy_upper_bound(p));
    printf("The Cauchy lower bound is %lg\n", polynomial_cauchy_lower_bound(p));
    printf("The Lagrange upper bound is %lg\n", polynomial_lagrange_upper_bound(p));
    printf("The Lagrange lower bound is %lg\n", polynomial_lagrange_lower_bound(p));
    printf("The Cauchy upper quota is %lg\n", polynomial_cauchy_upper_quota(p));
    printf("The Cauchy lower quota is %lg\n", polynomial_cauchy_lower_quota(p));
    printf("The Kojima's upper bound is %lg\n", polynomial_kojima_upper_bound(p));
    printf("The Kojima's lower bound is %lg\n", polynomial_kojima_lower_bound(p));
    double min, max;
    polynomial_root_bounds(p, &min, &max);
    printf("\nSo, the roots of the polynomial have absolute values between %lg and %lg\n\n", min, max);
    vector_dealloc(&p);
    return EXIT_SUCCESS;
}
