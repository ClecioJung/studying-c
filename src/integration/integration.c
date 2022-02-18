#include "../../lib/integration.h"

#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define PRECISION 0.005

typedef double (*Integration_Method_Fn)(const integration_fn fn, const double start, const double end, const size_t segments);

typedef struct {
    const char *const name;
    const Integration_Method_Fn method_fn;
} Integration_Algorithm;

double f(const double x) {
    return (x * x);
}

int test_algorithm(const Integration_Algorithm algorithm) {
    const double start = 0.0;
    const double end = 1.0;
    const size_t num_of_segments = 10;
    const double expected_result = 1.0 / 3.0;
    const double result = algorithm.method_fn(f, start, end, num_of_segments);
    if (are_close(result, expected_result, PRECISION)) {
        printf("Method %s estimated the integral from %lg to %lg to be: %lg\n", algorithm.name, start, end, result);
        return EXIT_SUCCESS;
    } else {
        fprintf(stderr, "The method %s couldn't estimate correctly the provided integral!\n", algorithm.name);
        return EXIT_FAILURE;
    }
}

int main(void) {
    const Integration_Algorithm algorithms[] = {
        {"trapezoidal rule", trapezoidal_rule},
        {"Simpson rule", simpson_rule},
        {"Simpson rule 3/8", simpson_rule_38},
    };
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(algorithms[0]); i++) {
        if (test_algorithm(algorithms[i]) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
