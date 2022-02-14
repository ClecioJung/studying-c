// MIT License

// Copyright (c) 2022 CLECIO JUNG <clecio.jung@gmail.com>

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

//------------------------------------------------------------------------------
// HEADER
//------------------------------------------------------------------------------

#ifndef __ODE_H
#define __ODE_H

// Header file with useful functions for dealing with
// Ordinary Differential Equations

#include <stdlib.h>

typedef double (*ode_fn)(const double x, const double y);
typedef void (*ode_method)(const ode_fn fn, double *const x, double *const y, const double step);

void euler_method(const ode_fn fn, double *const x, double *const y, const double step);
void heun_method(const ode_fn fn, double *const x, double *const y, const double step);
void runge_kutta_3_method(const ode_fn fn, double *const x, double *const y, const double step);
void runge_kutta_4_method(const ode_fn fn, double *const x, double *const y, const double step);
void butcher_method(const ode_fn fn, double *const x, double *const y, const double step);
void advance_ode_by_steps(const ode_method method, const ode_fn fn, double *const x, double *const y, const double step, const size_t number_of_steps);
double ode_final_value(const ode_method method, const ode_fn fn, const double initial_x, const double final_x, const double initial_y, const size_t number_of_steps);

#endif  // __ODE_H

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------

#ifdef ODE_IMPLEMENTATION

void euler_method(const ode_fn fn, double *const x, double *const y, const double step) {
    if ((x != NULL) && (y != NULL)) {
        *y = *y + step * fn(*x, *y);
        *x = *x + step;
    }
}

void heun_method(const ode_fn fn, double *const x, double *const y, const double step) {
    if ((x != NULL) && (y != NULL)) {
        const double k = *y + step * fn(*x, *y);
        *y = *y + (fn(*x, *y) + fn(*x + step, k)) * (step / 2.0);
        *x = *x + step;
    }
}

// Third-Order Runge-Kutta Method
void runge_kutta_3_method(const ode_fn fn, double *const x, double *const y, const double step) {
    if ((x != NULL) && (y != NULL)) {
        const double k1 = step * fn(*x, *y);
        const double k2 = step * fn(*x + step / 2.0, *y + k1 / 2.0);
        const double k3 = step * fn(*x + step, *y - k1 + 2.0 * k2);
        *y = *y + (k1 + 4.0 * k2 + k3) / 6.0;
        *x = *x + step;
    }
}

// Fourth-Order Runge-Kutta Method
void runge_kutta_4_method(const ode_fn fn, double *const x, double *const y, const double step) {
    if ((x != NULL) && (y != NULL)) {
        const double k1 = step * fn(*x, *y);
        const double k2 = step * fn(*x + step / 2.0, *y + k1 / 2.0);
        const double k3 = step * fn(*x + step / 2.0, *y + k2 / 2.0);
        const double k4 = step * fn(*x + step, *y + k3);
        *y = *y + (k1 + 2.0 * k2 + 2.0 * k3 + k4) / 6.0;
        *x = *x + step;
    }
}

// Butcher’s (1964) fifth-order RK method
void butcher_method(const ode_fn fn, double *const x, double *const y, const double step) {
    if ((x != NULL) && (y != NULL)) {
        const double k1 = step * fn(*x, *y);
        const double k2 = step * fn(*x + step / 4.0, *y + k1 / 4.0);
        const double k3 = step * fn(*x + step / 4.0, *y + k1 / 8.0 + k2 / 8.0);
        const double k4 = step * fn(*x + step / 2.0, *y - k2 / 2.0 + k3);
        const double k5 = step * fn(*x + 3.0 / 4.0 * step, *y + 3.0 / 16.0 * k1 + 9.0 / 16.0 * k4);
        const double k6 = step * fn(*x + step, *y - 3.0 / 7.0 * k1 + 2.0 / 7.0 * k2 + 12.0 / 7.0 * k3 - 12.0 / 7.0 * k4 + 8.0 / 7.0 * k5);
        *y = *y + (7.0 * k1 + 32.0 * k3 + 12.0 * k4 + 32.0 * k5 + 7.0 * k6) / 90.0;
        *x = *x + step;
    }
}

void advance_ode_by_steps(const ode_method method, const ode_fn fn, double *const x, double *const y, const double step, const size_t number_of_steps) {
    for (size_t i = 0; i < number_of_steps; i++) {
        method(fn, x, y, step);
    }
}

double ode_final_value(const ode_method method, const ode_fn fn, const double initial_x, const double final_x, const double initial_y, const size_t number_of_steps) {
    double x = initial_x;
    double y = initial_y;
    const double step = (final_x - initial_x) / ((double)number_of_steps);
    advance_ode_by_steps(method, fn, &x, &y, step, number_of_steps);
    return y;
}

#endif  // ODE_IMPLEMENTATION

//------------------------------------------------------------------------------
// END
//------------------------------------------------------------------------------

// MIT License

// Copyright (c) 2022 CLECIO JUNG <clecio.jung@gmail.com>

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.