#include <stdlib.h>
#include <stdio.h>

#define ODE_IMPLEMENTATION
#include "ode.h"

double f(const double x, const double y) {
    (void)x;
    return -y;
}

int main(void) {
    const size_t number_of_steps = 10;
    const double initial_x = 0.0;
    const double final_x = 1.0;
    const double initial_y = 1.0;
    printf("The final value of the ODE is: %g (Euler method)\n", ode_final_value(euler_method, f, initial_x, final_x, initial_y, number_of_steps));
    printf("The final value of the ODE is: %g (Heun method)\n", ode_final_value(heun_method, f, initial_x, final_x, initial_y, number_of_steps));
    printf("The final value of the ODE is: %g (Third-Order Runge-Kutta method)\n", ode_final_value(runge_kutta_3_method, f, initial_x, final_x, initial_y, number_of_steps));
    printf("The final value of the ODE is: %g (Fourth-Order Runge-Kutta method)\n", ode_final_value(runge_kutta_4_method, f, initial_x, final_x, initial_y, number_of_steps));
    printf("The final value of the ODE is: %g (Butcher's method)\n", ode_final_value(butcher_method, f, initial_x, final_x, initial_y, number_of_steps));
    return EXIT_SUCCESS;
}
