#include "../../lib/ode.h"

#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define PRECISION 0.02

typedef void (*ODE_Method_Fn)(const ode_fn fn, double *const x, double *const y, const double step);

typedef struct {
    const char *const name;
    const ODE_Method_Fn method_fn;
} ODE_Algorithm;

double f(const double x, const double y) {
    (void)x;
    return -y;
}

int test_algorithm(const ODE_Algorithm algorithm) {
    const size_t number_of_steps = 10;
    const double initial_x = 0.0;
    const double final_x = 1.0;
    const double initial_y = 1.0;
    const double expected_result = 0.367;
    const double result = ode_final_value(algorithm.method_fn, f, initial_x, final_x, initial_y, number_of_steps);
    if (are_close(result, expected_result, PRECISION)) {
        printf("Method %s estimated the final value of the ODE to be: %lg\n", algorithm.name, result);
        return EXIT_SUCCESS;
    } else {
        fprintf(stderr, "The method %s couldn't estimate correctly the provided ODE!\n", algorithm.name);
        return EXIT_FAILURE;
    }
}

int main(void) {
    const ODE_Algorithm algorithms[] = {
        {"Euler method", euler_method},
        {"Heun method", heun_method},
        {"Runge Kutta third method", runge_kutta_3_method},
        {"Runge Kutta fourth method", runge_kutta_4_method},
        {"Butcher method", butcher_method},
    };
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(algorithms[0]); i++) {
        if (test_algorithm(algorithms[i]) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}
