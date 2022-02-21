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

#include "curve-fitting.h"

#include <stdlib.h>

#include "linear-systems.h"
#include "matrix.h"
#include "scalar.h"

double lagrange_interpolation(const Vector x, const Vector y, const double value) {
    if (x.len != y.len) {
        return 0.0;  // Invalid operation
    }
    double sum = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        double product = y.data[i];
        for (size_t j = 0; j < x.len; j++) {
            if (j != i) {
                product *= (value - x.data[j]) / (x.data[i] - x.data[j]);
            }
        }
        sum += product;
    }
    return sum;
}

double linear_regression(const Vector x, const Vector y, double *const a, double *const b) {
    if ((x.len != y.len) || (a == NULL) || (b == NULL)) {
        return 0.0;  // Invalid operation
    }
    double sum_x = 0.0, sum_x_squared = 0.0;
    double sum_y = 0.0, sum_y_squared = 0.0;
    double sum_xy = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        sum_x += x.data[i];
        sum_x_squared += x.data[i] * x.data[i];
        sum_y += y.data[i];
        sum_y_squared += y.data[i] * y.data[i];
        sum_xy += x.data[i] * y.data[i];
    }
    *a = (x.len * sum_xy - sum_x * sum_y) / (x.len * sum_x_squared - sum_x * sum_x);
    *b = (sum_y * sum_x_squared - sum_xy * sum_x) / (x.len * sum_x_squared - sum_x * sum_x);
    double r = (x.len * sum_xy - sum_x * sum_y) / (square_root(x.len * sum_x_squared - sum_x * sum_x) * square_root(x.len * sum_y_squared - sum_y * sum_y));
    return r;
}

// Remember to free the returned vector after calling this function!
Vector polynomial_regression(const Vector x, const Vector y, const size_t order) {
    if ((x.len != y.len) || (x.len < (order + 1))) {
        return (Vector){0};  // Invalid operation
    }
    Matrix A = matrix_init(order + 1, order + 1, 0.0);
    Vector b = vector_init(order + 1, 0.0);
    if (!matrix_is_valid(A) || !vector_is_valid(b)) {
        matrix_dealloc(&A);
        vector_dealloc(&b);
        return (Vector){0};
    }
    // Computing augmented matrix
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j <= i; j++) {
            for (size_t l = 0; l < x.len; l++) {
                matrix_inc(A, i, j, power(x.data[l], (i + j)));
            }
            if (i != j) {
                matrix_set(A, j, i, matrix_get(A, i, j));
            }
        }
        for (size_t l = 0; l < x.len; l++) {
            b.data[i] += y.data[l] * power(x.data[l], i);
        }
    }
    // If the system doesn't have a siingle unique solution,
    // the coefficients vector will be deallocated automatically by the function gaussian_elimination()
    Vector coefficients;
    gaussian_elimination(A, b, &coefficients);
    matrix_dealloc(&A);
    vector_dealloc(&b);
    return coefficients;
}

double compute_polynomial(const Vector coefficients, const double x) {
    double y = 0.0;
    for (size_t i = 0; i < coefficients.len; i++) {
        y += coefficients.data[i] * power(x, i);
    }
    return y;
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