#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "../../lib/curve-fitting.h"
#include "../../lib/polynomial.h"
#include "../../lib/scalar.h"
#include "../../lib/vector.h"

#define PRECISION 1e-8

typedef Vector (*Polynomial_Interpolation)(const Vector, const Vector, const Vector);

typedef struct {
    const char *const name;
    const Polynomial_Interpolation interpolation;
} Compute_Algorithm;

int test_algorithm(const Compute_Algorithm algorithm, const Vector x, const Vector y) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    Vector result = algorithm.interpolation(x, y, x);
    gettimeofday(&stop, NULL);
    unsigned long int delta_us = (stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec);
    if (!vector_is_valid(result)) {
        fprintf(stderr, "Algorithm %s failed to allocate memory.\n", algorithm.name);
        return EXIT_FAILURE;
    }
    for (size_t i = 0; i < x.len; i++) {
        if (!are_close(result.data[i], y.data[i], PRECISION)) {
            fprintf(stderr, "Algorithm %s failed to correctly interpolate the given points.\n", algorithm.name);
            return EXIT_FAILURE;
        }
    }
    printf("Took %04ld us for the %s to compute %ld points.\n", delta_us, algorithm.name, x.len);
    vector_dealloc(&result);
    return EXIT_SUCCESS;
}

int main(void) {
    const Compute_Algorithm algorithms[] = {
        {"Polynomial interpolation", polynomial_interpolation_vector},
        {"Lagrange interpolation", lagrange_interpolation_vector},
        {"Gregory-Newton interpolation", gregory_newton_interpolation_vector},
    };
    const size_t len = 20;
    Vector x = vector_alloc(len);
    Vector y = vector_alloc(len);
    {
        Vector polynomial = vector_alloc(len);
        for (size_t i = 0; i < len; i++) {
            polynomial.data[i] = 1.0 - ((double)i) / ((double)len);
        }
        for (size_t i = 0; i < len; i++) {
            x.data[i] = ((double)i) / ((double)len);
            y.data[i] = polynomial_horner_evaluation(polynomial, x.data[i]);
        }
        vector_dealloc(&polynomial);
    }
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(algorithms[0]); i++) {
        if (test_algorithm(algorithms[i], x, y) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    vector_dealloc(&x);
    vector_dealloc(&y);
    return EXIT_SUCCESS;
}
