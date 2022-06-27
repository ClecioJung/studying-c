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
#include <stdio.h>
#include <stdlib.h>

#include "linear-systems.h"
#include "matrix.h"
#include "polynomial.h"
#include "scalar.h"
#include "vector.h"

#define PRECISION 1e-12

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

// Remember to free the returned spline after calling this function!
Spline spline_alloc(const size_t len) {
    Spline spline = (Spline){0};
    if (len != 0) {
        spline.x = (double *)malloc((len + 1) * sizeof(double));
        spline.coefficients = (double(*)[SPLINE_ORDER + 1]) malloc(len * sizeof(*spline.coefficients));
        if ((spline.x == NULL) || (spline.coefficients == NULL)) {
            free(spline.x);
            free(spline.coefficients);
            spline.x = NULL;
            spline.coefficients = NULL;
        } else {
            spline.len = len;
        }
    }
    return spline;
}

void spline_dealloc(Spline *const spline) {
    if (spline == NULL) {
        return;
    }
    free(spline->x);
    free(spline->coefficients);
    spline->x = NULL;
    spline->coefficients = NULL;
    spline->len = 0;
}

// Remember to free the returned spline after calling this function!
// It is expected that the data vector x is in ascending order!
// S_start and S_end are only used if spline_type == Spline_Clamped
Spline spline_interpolation(const Vector x, const Vector y, const Spline_Type spline_type, const double S_start, const double S_end) {
    if ((x.len <= 1) || (x.len != y.len)) {
        return (Spline){0};  // Invalid operation
    }
    Spline spline = spline_alloc(x.len - 1);
    Vector t = vector_alloc(x.len);
    Vector r = vector_alloc(x.len);
    Vector d = vector_alloc(x.len);
    Vector S = vector_alloc(x.len);
    if ((spline.len == 0) || !vector_is_valid(t) || !vector_is_valid(r) || !vector_is_valid(d) || !vector_is_valid(S)) {
        spline_dealloc(&spline);
        vector_dealloc(&t);
        vector_dealloc(&r);
        vector_dealloc(&d);
        vector_dealloc(&S);
        return (Spline){0};
    }
    {  // Coefficients not used in tridiagonal systems solving algorithm
        t.data[0] = 0.0;
        d.data[x.len - 1] = 0.0;
    }
    double previous_h = x.data[1] - x.data[0];
    for (size_t i = 1; (i + 1) < x.len; i++) {
        const double h = x.data[i + 1] - x.data[i];
        t.data[i] = previous_h;
        r.data[i] = 2.0 * (previous_h + h);
        d.data[i] = h;
        S.data[i] = 6.0 * ((y.data[i + 1] - y.data[i]) / h - (y.data[i] - y.data[i - 1]) / previous_h);
        previous_h = h;
    }
    switch (spline_type) {
        // Clamped spline => S_start and S_end known
        case Spline_Clamped:
            r.data[0] = 1.0;
            d.data[0] = 0.0;
            S.data[0] = S_start;
            t.data[x.len - 1] = 0.0;
            r.data[x.len - 1] = 1.0;
            S.data[x.len - 1] = S_end;
            break;
        // Quadratic spline => S[0] = S[1] and S[end] = S[end-1]
        case Spline_Quadratic:
            r.data[0] = 1.0;
            d.data[0] = -1.0;
            S.data[0] = 0.0;
            t.data[x.len - 1] = -1.0;
            r.data[x.len - 1] = 1.0;
            S.data[x.len - 1] = 0.0;
            break;
        // Natural spline => S_start = S_end = 0.0
        case Spline_Natural:
        default:
            r.data[0] = 1.0;
            d.data[0] = 0.0;
            S.data[0] = 0.0;
            t.data[x.len - 1] = 0.0;
            r.data[x.len - 1] = 1.0;
            S.data[x.len - 1] = 0.0;
            break;
    }
    tridiagonal_solving_over(t, r, d, S);
    for (size_t i = 0; (i + 1) < x.len; i++) {
        spline.x[i] = x.data[i];
        const double h = x.data[i + 1] - x.data[i];
        // coefficient a
        spline.coefficients[i][3] = (S.data[i + 1] - S.data[i]) / (6.0 * h);
        // coefficient b
        spline.coefficients[i][2] = S.data[i] / 2.0;
        // coefficient c
        spline.coefficients[i][1] = (y.data[i + 1] - y.data[i]) / h - (S.data[i + 1] + 2.0 * S.data[i]) * h / 6.0;
        // coefficient d
        spline.coefficients[i][0] = y.data[i];
    }
    spline.x[x.len - 1] = x.data[x.len - 1];
    vector_dealloc(&t);
    vector_dealloc(&r);
    vector_dealloc(&d);
    vector_dealloc(&S);
    return spline;
}

void spline_print(const Spline spline) {
    for (size_t i = 0; i < spline.len; i++) {
        const double x = spline.x[i];
        printf("[%03zu] x in (%lg - %lg):", i, x, spline.x[i + 1]);
        size_t printed_terms = 0;
        for (size_t j = SPLINE_ORDER; j <= SPLINE_ORDER; j--) {
            const double coef = spline.coefficients[i][j];
            if (are_close(coef, 0.0, PRECISION)) {
                continue;
            }
            if (coef < 0) {
                printf(" - ");
            } else if (printed_terms > 0) {
                printf(" + ");
            } else {
                printf(" ");
            }
            printf("%lg", fabs(coef));
            if (j != 0) {
                if (!are_close(coef, 0.0, PRECISION)) {
                    printf("*");
                }
                if (!are_close(x, 0.0, PRECISION)) {
                    printf("(x %c %lg)", ((x < 0) ? '+' : '-'), fabs(x));
                } else {
                    printf("x");
                }
                if (j != 1) {
                    printf("^%ld", j);
                }
                printed_terms++;
            }
        }
        if (printed_terms == 0) {
            printf("0");
        }
        printf("\n");
    }
}

size_t spline_range_search(const Spline spline, const double value) {
    // Value is otside of the spline range
    if ((value < spline.x[0]) || (value > spline.x[spline.len])) {
        return spline.len;
    }
    // Use binary search to locate the value in the ranges of the vector x
    size_t middle = spline.len;
    size_t low = 0;
    size_t high = spline.len;
    while (low <= high) {
        middle = (low + high) / 2;
        if (value < spline.x[middle]) {
            high = middle - 1;
        } else if (value > spline.x[middle]) {
            low = middle + 1;
        } else {
            break;
        }
    }
    if ((middle != 0) && (value > spline.x[middle - 1])) {
        return (middle - 1);
    }
    return middle;
}

double spline_evaluation(const Spline spline, const double value) {
    const size_t index = spline_range_search(spline, value);
    if (index >= spline.len) {
        return NAN;  // Invalid range
    }
    // Polynomial evaluation using Horner's method
    double y = 0.0;
    const double x = value - spline.x[index];
    for (size_t j = SPLINE_ORDER; j <= SPLINE_ORDER; j--) {
        y = spline.coefficients[index][j] + y * x;
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