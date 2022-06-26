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

#include <math.h>
#include <stdlib.h>

#include "linear-systems.h"
#include "matrix.h"
#include "polynomial.h"
#include "scalar.h"
#include "vector.h"

double polynomial_interpolation(const Vector x, const Vector y, const double value) {
    if (x.len != y.len) {
        return NAN;  // Invalid operation
    }
    Matrix A = matrix_vandermonde(x, x.len);
    Vector b = vector_copy(y);
    if (!matrix_is_valid(A) || !vector_is_valid(b)) {
        matrix_dealloc(&A);
        vector_dealloc(&b);
        return NAN;
    }
    gaussian_elimination_over(A, b);
    const double result = polynomial_horner_evaluation(b, value);
    matrix_dealloc(&A);
    vector_dealloc(&b);
    return result;
}

double lagrange_interpolation(const Vector x, const Vector y, const double value) {
    if (x.len != y.len) {
        return NAN;  // Invalid operation
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

// Remember to free the returned vector after calling this function!
Vector gregory_newton_differences(const Vector x, const Vector y) {
    if (x.len != y.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector differences = vector_alloc(x.len);
    Vector diff_auxiliar = vector_copy(y);
    if (!vector_is_valid(differences) || !vector_is_valid(diff_auxiliar)) {
        vector_dealloc(&differences);
        vector_dealloc(&diff_auxiliar);
        return (Vector){0};
    }
    for (size_t i = 0; (i + 1) < x.len; i++) {
        for (size_t j = 0; (j + i + 1) < x.len; j++) {
            diff_auxiliar.data[j] = (diff_auxiliar.data[j + 1] - diff_auxiliar.data[j]) / (x.data[j + i + 1] - x.data[j]);
        }
        differences.data[i] = diff_auxiliar.data[0];
    }
    differences.data[x.len - 1] = 0.0;
    vector_dealloc(&diff_auxiliar);
    return differences;
}

double gregory_newton_interpolation(const Vector x, const Vector y, const double value) {
    if (x.len != y.len) {
        return NAN;  // Invalid operation
    }
    Vector differences = gregory_newton_differences(x, y);
    if (!vector_is_valid(differences)) {
        return NAN;
    }
    double sum = y.data[0];
    double product = 1.0;
    for (size_t i = 0; i < x.len; i++) {
        product *= value - x.data[i];
        sum += differences.data[i] * product;
    }
    vector_dealloc(&differences);
    return sum;
}

double linear_regression(const Vector x, const Vector y, double *const a, double *const b) {
    if ((x.len != y.len) || (a == NULL) || (b == NULL)) {
        return NAN;  // Invalid operation
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
    const double r2 = power((x.len * sum_xy - sum_x * sum_y), 2) / ((x.len * sum_x_squared - sum_x * sum_x) * (x.len * sum_y_squared - sum_y * sum_y));
    return r2;
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
    gaussian_elimination_over(A, b);
    matrix_dealloc(&A);
    return b;
}

double r_squared(const Vector x, const Vector y, const Vector polynomial) {
    if (x.len != y.len) {
        return NAN;  // Invalid operation
    }
    const double y_mean = vector_mean(y);
    double ssr = 0.0;  // Residual sum of squares
    double sst = 0.0;  // Total sum of squares
    for (size_t i = 0; i < x.len; i++) {
        const double y_estimated = polynomial_horner_evaluation(polynomial, x.data[i]);
        ssr += power((y.data[i] - y_estimated), 2);
        sst += power((y.data[i] - y_mean), 2);
    }
    return (1.0 - (ssr / sst));
}

// Remember to free the returned vector after calling this function!
Vector polynomial_interpolation_vector(const Vector x, const Vector y, const Vector values) {
    if (x.len != y.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector result = vector_init(values.len, 0.0);
    Matrix A = matrix_vandermonde(x, x.len);
    Vector b = vector_copy(y);
    if (!vector_is_valid(result) || !matrix_is_valid(A) || !vector_is_valid(b)) {
        vector_dealloc(&result);
        matrix_dealloc(&A);
        vector_dealloc(&b);
        return (Vector){0};
    }
    lu_solving_over(A, b);
    matrix_dealloc(&A);
    for (size_t i = 0; i < values.len; i++) {
        result.data[i] = polynomial_horner_evaluation(b, values.data[i]);
    }
    vector_dealloc(&b);
    return result;
}

// Remember to free the returned vector after calling this function!
Vector lagrange_interpolation_vector(const Vector x, const Vector y, const Vector values) {
    if (x.len != y.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector result = vector_init(values.len, 0.0);
    if (!vector_is_valid(result)) {
        return (Vector){0};
    }
    for (size_t i = 0; i < x.len; i++) {
        double product = y.data[i];
        for (size_t j = 0; j < x.len; j++) {
            if (j != i) {
                product *= (values.data[i] - x.data[j]) / (x.data[i] - x.data[j]);
            }
        }
        result.data[i] += product;
    }
    return result;
}

Vector gregory_newton_interpolation_vector(const Vector x, const Vector y, const Vector values) {
    if (x.len != y.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector result = vector_init(values.len, y.data[0]);
    Vector differences = gregory_newton_differences(x, y);
    if (!vector_is_valid(result) || !vector_is_valid(differences)) {
        vector_dealloc(&result);
        vector_dealloc(&differences);
        return (Vector){0};
    }
    for (size_t i = 0; i < x.len; i++) {
        double product = 1.0;
        for (size_t j = 0; j < x.len; j++) {
            product *= values.data[i] - x.data[j];
            result.data[i] += differences.data[j] * product;
        }
    }
    vector_dealloc(&differences);
    return result;
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