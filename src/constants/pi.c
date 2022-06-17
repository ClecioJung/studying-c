#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "../../lib/scalar.h"

#define MAX_ITERATIONS 1000000000
#define PRECISION 1e-8

double pi_Gregory_leibniz(void) {
    double pi = 0.0;
    for (uint64_t i = 0; i < MAX_ITERATIONS; i++) {
        // Gregory-Leibniz Series
        const double delta = (i % 2 ? -4.0 : 4.0) / ((double)(2 * i + 1));
        pi += delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return pi;
}

double pi_Nilakantha(void) {
    double pi = 3.0;
    for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
        // Nilakantha Series
        const double delta = (i % 2 ? 4.0 : -4.0) / ((double)((2 * i) * (2 * i + 1) * (2 * i + 2)));
        pi += delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return pi;
}

double pi_Ramanujan(void) {
    double pi = 0.0;
    for (uint64_t i = 0; i < MAX_ITERATIONS; i++) {
        const double previous_pi = pi;
        // Srinivasa Ramanujan formula based on modular equations
        const double inverse_delta = (2.0 * sqrt(2.0) / 9801.0) * ((double)(factorial(4 * i) * (1103 + 26390 * i))) /
                                     ((double)(power(factorial(i), 4) * power(396, 4 * i)));
        pi = (previous_pi > 0) ? (1.0 / (1.0 / previous_pi + inverse_delta)) : (1.0 / inverse_delta);
        if (fabs(pi - previous_pi) < PRECISION) {
            break;
        }
    }
    return pi;
}

double pi_Chudnovsky(void) {
    double pi = 0.0;
    for (uint64_t i = 0; i < MAX_ITERATIONS; i++) {
        const double previous_pi = pi;
        // Chudnovsky formula
        const double inverse_delta = (12.0 / pow(640320.0, 1.5)) * ((double)(factorial(6 * i) * (13591409 + 545140134 * i))) /
                                     ((double)(factorial(3 * i) * power(factorial(i), 3) * power(-640320, 3 * i)));
        pi = (previous_pi > 0) ? (1.0 / (1.0 / previous_pi + inverse_delta)) : (1.0 / inverse_delta);
        if (fabs(pi - previous_pi) < PRECISION) {
            break;
        }
    }
    return pi;
}

typedef double (*Approx_Fn)(void);

typedef struct {
    const char *const name;
    const Approx_Fn approx_fn;
} Approx_Algorithm;

int test_algorithm(const Approx_Algorithm algorithm) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    const double result = algorithm.approx_fn();
    gettimeofday(&stop, NULL);
    unsigned long int delta_us = (stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec);
    if (are_close(result, PI, PRECISION)) {
        printf("Took %04ld us for the %s to approximate pi to be %lg\n", delta_us, algorithm.name, result);
        return EXIT_SUCCESS;
    } else {
        fprintf(stderr, "Algorithm %s failed to approximate pi: %lg\n", algorithm.name, result);
        return EXIT_FAILURE;
    }
}

int main(void) {
    const Approx_Algorithm algorithms[] = {
        {"Gregory-Leibniz Series", pi_Gregory_leibniz},
        {"Nilakantha Series", pi_Nilakantha},
        {"Srinivasa Ramanujan formula", pi_Ramanujan},
        {"Chudnovsky formula", pi_Chudnovsky},
    };
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(algorithms[0]); i++) {
        if (test_algorithm(algorithms[i]) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
