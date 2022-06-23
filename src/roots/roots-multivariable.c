#include <stdio.h>
#include <stdlib.h>

#include "../../lib/roots.h"
#include "../../lib/scalar.h"
#include "../../lib/vector.h"

// Remember to free the returned vector after calling this function!
Vector f(const Vector x) {
    if (x.len != 2) {
        return (Vector){0};  // Invalid operation
    }
    Vector y = vector_alloc(2);
    if (vector_is_valid(y)) {
        y.data[0] = exponential(x.data[0]) + x.data[1] - 1.0;
        y.data[1] = x.data[0] * x.data[0] + x.data[1] * x.data[1] - 4.0;
    }
    return y;
}

int main(void) {
    Vector initial = vector_alloc(2);
    initial.data[0] = 1.0;
    initial.data[1] = -1.0;
    {
        Vector roots = multivariable_newton_method(f, initial);
        Vector image = f(roots);
        if (vector_is_null(image)) {
            printf("Roots found by the Newton-Raphson method:\n");
            vector_print(roots);
        } else {
            fprintf(stderr, "The Newton-Raphson method didn't found the correct roots!\n");
            return EXIT_FAILURE;
        }
        vector_dealloc(&roots);
        vector_dealloc(&image);
    }
    {
        Vector roots = multivariable_broyden_method(f, initial);
        Vector image = f(roots);
        if (vector_is_null(image)) {
            printf("Roots found by the Broyden's method:\n");
            vector_print(roots);
        } else {
            fprintf(stderr, "The Broyden's method didn't found the correct roots!\n");
            return EXIT_FAILURE;
        }
        vector_dealloc(&roots);
        vector_dealloc(&image);
    }
    vector_dealloc(&initial);
    return EXIT_SUCCESS;
}
