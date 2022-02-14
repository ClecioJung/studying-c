#include <stdio.h>
#include <stdlib.h>

#define ROOTS_IMPLEMENTATION
#include "roots.h"

double f(const double x) {
    return (x * x - 2.0);
}

double df(const double x) {
    return (2.0 * x);
}

int main(void) {
    printf("Root found by the bissection method: %g\n", bisection_method(f, 0.0, 2.0));
    printf("Root found by the fake-position method: %g\n", fakepos_method(f, 0.0, 2.0));
    printf("Root found by the Newton-Raphson method: %g\n", newton_raphson_method(f, df, 2.0));
    printf("Root found by the secant method: %g\n", secant_method(f, 0.0, 2.0));
    return EXIT_SUCCESS;
}
