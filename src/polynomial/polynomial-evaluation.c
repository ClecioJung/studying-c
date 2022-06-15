#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "../../lib/polynomial.h"
#include "../../lib/scalar.h"
#include "../../lib/vector.h"

#define PRECISION 1e-8

typedef double (*Compute_Pol)(const Vector, const double);

typedef struct {
    const char *const name;
    const Compute_Pol compute_pol;
} Compute_Algorithm;

int test_algorithm(const Compute_Algorithm algorithm, const Vector coefficients, const double x, const double expected_result) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    const double result = algorithm.compute_pol(coefficients, x);
    gettimeofday(&stop, NULL);
    unsigned long int delta_us = (stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec);
    if (are_close(result, expected_result, PRECISION)) {
        printf("Took %04ld us for the %s to compute the given polynomial at %lg to be: %lg\n", delta_us, algorithm.name, x, result);
        return EXIT_SUCCESS;
    } else {
        fprintf(stderr, "Algorithm %s failed to compute the given polynomial: %lg\n", algorithm.name, result);
        return EXIT_FAILURE;
    }
}

int main(void) {
    const Compute_Algorithm algorithms[] = {
        {"conventional method", polynomial_evaluation},
        {"Horner's method", polynomial_horner_evaluation},
        {"Ruffini's method", polynomial_ruffini_evaluation},
    };
    const size_t len = 10000;
    const double coefficients_value = 1.0;
    const double x = 1.0;
    const double y = (double)len;
    const Vector coefficients = vector_init(len, coefficients_value);
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(algorithms[0]); i++) {
        if (test_algorithm(algorithms[i], coefficients, x, y) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
