#include <stdio.h>
#include <stdlib.h>

#include "../../lib/curve-fitting.h"
#include "../../lib/scalar.h"

#define PRECISION 1e-10

int main(void) {
    Vector x = vector_alloc(6);
    Vector y = vector_alloc(6);
    x.data[0] = 0.0;
    x.data[1] = 1.0;
    x.data[2] = 2.0;
    x.data[3] = 3.0;
    x.data[4] = 4.0;
    x.data[5] = 5.0;
    y.data[0] = 0.0;
    y.data[1] = 2.0;
    y.data[2] = 4.0;
    y.data[3] = 6.0;
    y.data[4] = 2.0;
    y.data[5] = 0.0;
    printf("Polynomial interpolation:\n");
    for (size_t i = 0; i < x.len; i++) {
        const double result = polynomial_interpolation(x, y, x.data[i]);
        if (are_close(result, y.data[i], PRECISION)) {
            printf("Interpolation at %lg results in %lg\n", x.data[i], result);
        } else {
            fprintf(stderr, "Polynomial interpolation: couldn't interpolate correctly the provided values!\n");
            return EXIT_FAILURE;
        }
    }
    vector_dealloc(&x);
    vector_dealloc(&y);
    return EXIT_SUCCESS;
}
