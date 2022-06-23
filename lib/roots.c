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
// SOURCE
//------------------------------------------------------------------------------

#include "roots.h"

#include <math.h>
#include <stdlib.h>

#include "linear-systems.h"
#include "matrix.h"
#include "scalar.h"
#include "vector.h"

#define MAX_ITERATIONS 10000
#define PRECISION 1e-10

double bisection_method(const root_fn fn, double start, double end) {
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double middle = (start + end) / 2.0;
        if (fabs(end - start) < PRECISION) {
            return middle;
        }
        if (fn(start) * fn(middle) < 0) {
            end = middle;
        } else {
            start = middle;
        }
    }
    return 0.0;
}

double fakepos_method(const root_fn fn, double start, double end) {
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double x = end - ((fn(end) * (end - start)) / (fn(end) - fn(start)));
        if (fabs(end - start) < PRECISION) {
            return x;
        }
        if (fn(start) * fn(end) < 0) {
            end = x;
        } else {
            start = x;
        }
    }
    return 0.0;
}

double newton_raphson_method(const root_fn fn, const root_fn dfn, const double initial) {
    double x = initial;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double delta = fn(x) / dfn(x);
        x -= delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return x;
}

double secant_method(const root_fn fn, const double start, const double end) {
    double previous_x = start;
    double x = end;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double delta = (fn(x) * (x - previous_x)) / (fn(x) - fn(previous_x));
        previous_x = x;
        x -= delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return x;
}

double muller_method(const root_fn fn, const double initial) {
    double x = initial;
    double previous_x = 0.99 * initial;
    double before_previous_x = 0.98 * initial;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double q = (x - previous_x) / (previous_x - before_previous_x);
        const double a = q * fn(x) - q * (1.0 + q) * fn(previous_x) + q * q * fn(before_previous_x);
        const double b = (2.0 * q + 1.0) * fn(x) - (1.0 + q) * (1.0 + q) * fn(previous_x) + q * q * fn(before_previous_x);
        const double c = (1.0 + q) * fn(x);
        const double delta = (2.0 * c * (x - previous_x)) / (b + sign(b) * square_root(b * b - 4.0 * a * c));
        before_previous_x = previous_x;
        previous_x = x;
        x -= delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return x;
}

// Remember to free the returned vector after calling this function!
Vector multivariable_newton_method(const multivariable_root_fn fn, const Vector initial) {
    Vector x = vector_copy(initial);
    Matrix J = matrix_alloc(x.len, x.len);
    if ((!vector_is_valid(x)) || (!matrix_is_valid(J))) {
        vector_dealloc(&x);
        matrix_dealloc(&J);
        return (Vector){0};
    }
    const double increment = PRECISION;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        Vector y = fn(x);
        if (!vector_is_valid(y)) {
            vector_dealloc(&x);
            matrix_dealloc(&J);
            return (Vector){0};
        }
        // Computes the Jacobian matrix J
        for (size_t j = 0; j < J.cols; j++) {
            // Increment x_j in order to calculate derivative
            x.data[j] += increment;
            Vector yj = fn(x);
            if (!vector_is_valid(yj)) {
                vector_dealloc(&y);
                vector_dealloc(&x);
                matrix_dealloc(&J);
                return (Vector){0};
            }
            for (size_t i = 0; i < J.rows; i++) {
                const double diff = (yj.data[i] - y.data[i]) / increment;
                matrix_set(J, i, j, diff);
            }
            x.data[j] -= increment;  // Undo increment
            vector_dealloc(&yj);
        }
        // Saves delta_x in y
        gaussian_elimination_over(J, y);
        vector_sub_over(x, y);
        const double delta_norm = vector_max_abs(y);
        vector_dealloc(&y);
        if (delta_norm < PRECISION) {
            break;
        }
    }
    matrix_dealloc(&J);
    return x;
}

// Remember to free the returned vector after calling this function!
Vector multivariable_broyden_method(const multivariable_root_fn fn, const Vector initial) {
    Vector x = vector_copy(initial);
    Matrix J = matrix_jacobian(fn, x);
    Matrix Jinv = matrix_inverse(J);
    matrix_dealloc(&J);
    if ((!vector_is_valid(x)) || (!matrix_is_valid(Jinv))) {
        vector_dealloc(&x);
        matrix_dealloc(&Jinv);
        return (Vector){0};
    }
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        Vector y = fn(x);
        Vector delta_x = matrix_mul_vector(Jinv, y);
        if ((!vector_is_valid(y)) || (!vector_is_valid(delta_x))) {
            vector_dealloc(&delta_x);
            vector_dealloc(&y);
            vector_dealloc(&x);
            matrix_dealloc(&Jinv);
            return (Vector){0};
        }
        vector_sub_over(x, delta_x);
        const double delta_norm = vector_max_abs(delta_x);
        {  // Uses Sherman–Morrison formula to compute the Jacobian inverse for the next iteration
            Vector J_times_y = matrix_mul_vector(Jinv, y);
            Vector dx_less_Jy = vector_sub(delta_x, J_times_y);
            if ((!vector_is_valid(J_times_y)) && (!vector_is_valid(dx_less_Jy))) {
                vector_dealloc(&dx_less_Jy);
                vector_dealloc(&J_times_y);
                vector_dealloc(&delta_x);
                vector_dealloc(&x);
                matrix_dealloc(&Jinv);
                return (Vector){0};
            }
            const double num = vector_dot_product(dx_less_Jy, delta_x);
            const double den = vector_dot_product(delta_x, J_times_y);
            matrix_scale_over((1.0 + num / den), Jinv);
            vector_dealloc(&J_times_y);
            vector_dealloc(&dx_less_Jy);
        }
        vector_dealloc(&y);
        vector_dealloc(&delta_x);
        if (delta_norm < PRECISION) {
            break;
        }
    }
    matrix_dealloc(&Jinv);
    return x;
}

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