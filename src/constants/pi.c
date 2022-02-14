#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define MATRIX_IMPLEMENTATION
#include "../math-vector.h"

data_type pi_Gregory_leibniz(void) {
    data_type pi = 0.0;
    for (uint64_t i = 0; i < MAX_ITERATIONS; i++) {
        // Gregory-Leibniz Series
        const data_type delta = (i % 2 ? -4.0 : 4.0) / ((data_type)(2 * i + 1));
        pi += delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return pi;
}

data_type pi_Nilakantha(void) {
    data_type pi = 3.0;
    for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
        // Nilakantha Series
        const data_type delta = (i % 2 ? 4.0 : -4.0) / ((data_type)((2 * i) * (2 * i + 1) * (2 * i + 2)));
        pi += delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return pi;
}

data_type pi_Ramanujan(void) {
    data_type pi = 0.0;
    for (uint64_t i = 0; i < MAX_ITERATIONS; i++) {
        const data_type previous_pi = pi;
        // Srinivasa Ramanujan formula based on modular equations
        const data_type inverse_delta = (2.0 * sqrt(2.0) / 9801.0) * ((data_type)(factorial(4 * i) * (1103 + 26390 * i))) /
                                     ((data_type)(power(factorial(i), 4) * power(396, 4 * i)));
        pi = (previous_pi > 0) ? (1.0 / (1.0 / previous_pi + inverse_delta)) : (1.0 / inverse_delta);
        if (fabs(pi - previous_pi) < PRECISION) {
            break;
        }
    }
    return pi;
}

data_type pi_Chudnovsky(void) {
    data_type pi = 0.0;
    for (uint64_t i = 0; i < MAX_ITERATIONS; i++) {
        const data_type previous_pi = pi;
        // Chudnovsky formula
        const data_type inverse_delta = (12.0 / pow(640320.0, 1.5)) * ((data_type)(factorial(6 * i) * (13591409 + 545140134 * i))) /
                                     ((data_type)(factorial(3 * i) * power(factorial(i), 3) * power(-640320, 3 * i)));
        pi = (previous_pi > 0) ? (1.0 / (1.0 / previous_pi + inverse_delta)) : (1.0 / inverse_delta);
        if (fabs(pi - previous_pi) < PRECISION) {
            break;
        }
    }
    return pi;
}

typedef data_type (*Approx_Fn)(void);

void test_algorithm(const char *const algorithm, const Approx_Fn approx_fn) {
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    const data_type result = approx_fn();
    gettimeofday(&stop, NULL);
    unsigned long int delta_us = (stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec);
    printf("Took %04ld us for the %s to approximate pi to be %lg\n", delta_us, algorithm, result);
}

int main(void) {
    test_algorithm("Gregory-Leibniz Series", pi_Gregory_leibniz);
    test_algorithm("Nilakantha Series", pi_Nilakantha);
    test_algorithm("Srinivasa Ramanujan formula", pi_Ramanujan);
    test_algorithm("Chudnovsky formula", pi_Chudnovsky);
    return EXIT_SUCCESS;
}
