#include "../../lib/roots.h"

#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define SQRT2 1.41421356237309504880
#define PRECISION 1e-10

double f(const double x) {
    return (x * x - 2.0);
}

double df(const double x) {
    return (2.0 * x);
}

int main(void) {
    {
        const double result = bisection_method(f, 0.0, 2.0);
        if (are_close(result, SQRT2, PRECISION)) {
            printf("Root found by the bissection method: %lg\n", result);
        } else {
            fprintf(stderr, "The bissection method didn't found the correct root\n");
            return EXIT_FAILURE;
        }
    }
    {
        const double result = fakepos_method(f, 0.0, 2.0);
        if (are_close(result, SQRT2, PRECISION)) {
            printf("Root found by the fake-position method: %lg\n", result);
        } else {
            fprintf(stderr, "The fake-position method didn't found the correct root\n");
            return EXIT_FAILURE;
        }
    }
    {
        const double result = newton_raphson_method(f, df, 2.0);
        if (are_close(result, SQRT2, PRECISION)) {
            printf("Root found by the Newton-Raphson method: %lg\n", result);
        } else {
            fprintf(stderr, "The Newton-Raphson method didn't found the correct root\n");
            return EXIT_FAILURE;
        }
    }
    {
        const double result = secant_method(f, 0.0, 2.0);
        if (are_close(result, SQRT2, PRECISION)) {
            printf("Root found by the secant method: %lg\n", result);
        } else {
            fprintf(stderr, "The secant method didn't found the correct root\n");
            return EXIT_FAILURE;
        }
    }
    {
        const double result = muller_method(f, 2.0);
        if (are_close(result, SQRT2, PRECISION)) {
            printf("Root found by the Muller method: %lg\n", result);
        } else {
            fprintf(stderr, "The Muller method didn't found the correct root\n");
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
