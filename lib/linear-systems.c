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

#include "linear-systems.h"

#include <math.h>
#include <stdbool.h>
#include <stdlib.h>

#include "scalar.h"

#define MAX_ITERATIONS 10000
#define PRECISION 1e-10

// Remember to free the returned vector after calling this function!
Vector back_substitution(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = vector_alloc(b.len);
    if (vector_is_valid(x)) {
        for (size_t i = (b.len - 1); i < b.len; i--) {
            double sum = b.data[i];
            for (size_t j = (i + 1); j < b.len; j++) {
                sum -= matrix_get(A, i, j) * x.data[j];
            }
            x.data[i] = sum / matrix_get(A, i, i);
        }
    }
    return x;
}

// Remember to free the returned vector after calling this function!
Vector forward_substitution(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = vector_alloc(b.len);
    if (vector_is_valid(x)) {
        for (size_t i = 0; i < b.len; i++) {
            double sum = b.data[i];
            for (size_t j = 0; j < i; j++) {
                sum -= matrix_get(A, i, j) * x.data[j];
            }
            x.data[i] = sum / matrix_get(A, i, i);
        }
    }
    return x;
}

// Remember to free the returned vector after calling this function!
Vector gaussian_elimination(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Matrix A_copied = matrix_copy(A);
    Vector b_copied = vector_copy(b);
    if (!matrix_is_valid(A_copied) || !vector_is_valid(b_copied)) {
        return (Vector){0};
    }
    // Partial pivoting
    for (size_t k = 0; (k + 1) < A_copied.rows; k++) {
        double w = fabs(matrix_get(A_copied, k, k));
        size_t r = k;
        for (size_t i = (k + 1); i < A_copied.rows; i++) {
            if (fabs(matrix_get(A_copied, i, k)) > w) {
                w = fabs(matrix_get(A_copied, i, k));
                r = i;
            }
        }
        if (r != k) {
            for (size_t i = k; i < A_copied.rows; i++) {
                const double temp = matrix_get(A_copied, k, i);
                matrix_set(A_copied, k, i, matrix_get(A_copied, r, i));
                matrix_set(A_copied, r, i, temp);
            }
            swap(&b_copied.data[k], &b_copied.data[r]);
        }
    }
    // Gaussian elimination
    for (size_t k = 0; (k + 1) < A_copied.rows; k++) {
        for (size_t i = (k + 1); i < A_copied.rows; i++) {
            const double m = matrix_get(A_copied, i, k) / matrix_get(A_copied, k, k);
            for (size_t j = (k + 1); j < A_copied.rows; j++) {
                matrix_dec(A_copied, i, j, m * matrix_get(A_copied, k, j));
            }
            b_copied.data[i] -= m * b_copied.data[k];
        }
    }
    Vector x = back_substitution(A_copied, b_copied);
    matrix_dealloc(&A_copied);
    vector_dealloc(&b_copied);
    return x;
}

// Remember to free the returned vector after calling this function!
Vector gauss_jordan(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Matrix A_copied = matrix_copy(A);
    Vector b_copied = vector_copy(b);
    if (!matrix_is_valid(A_copied) || !vector_is_valid(b_copied)) {
        return (Vector){0};
    }
    for (size_t k = 0; k < A_copied.rows; k++) {
        for (size_t i = 0; i < A_copied.rows; i++) {
            const double m = matrix_get(A_copied, i, k) / matrix_get(A_copied, k, k);
            const double p = matrix_get(A_copied, k, k);
            if (i == k) {
                for (size_t j = k; j < A_copied.cols; j++) {
                    matrix_set(A_copied, i, j, matrix_get(A_copied, i, j) / p);
                }
                b_copied.data[i] /= p;
            } else {
                for (size_t j = k; j < A_copied.cols; j++) {
                    matrix_dec(A_copied, i, j, m * matrix_get(A_copied, k, j));
                }
                b_copied.data[i] -= m * b_copied.data[k];
            }
        }
    }
    matrix_dealloc(&A_copied);
    return b_copied;
}

// Remember to free the returned vector after calling this function!
Vector lu_solving(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = (Vector){0};
    Matrix L, U;
    lu_decomposition(A, &L, &U);
    if (matrix_is_valid(L) || matrix_is_valid(U)) {
        Vector d = forward_substitution(L, b);
        if (vector_is_valid(d)) {
            x = back_substitution(U, d);
        }
        vector_dealloc(&d);
    }
    matrix_dealloc(&L);
    matrix_dealloc(&U);
    return x;
}

// Remember to free the returned vector after calling this function!
Vector jacobi_method(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = vector_init(b.len, 0.0);
    Vector previous_x = vector_alloc(b.len);
    if (!vector_is_valid(x) || !vector_is_valid(previous_x)) {
        vector_dealloc(&x);
    } else {
        for (size_t k = 0; k < MAX_ITERATIONS; k++) {
            vector_copy_over(&previous_x, x);
            for (size_t i = 0; i < x.len; i++) {
                double sum = b.data[i];
                for (size_t j = 0; j < x.len; j++) {
                    if (i != j) {
                        sum -= matrix_get(A, i, j) * previous_x.data[j];
                    }
                }
                x.data[i] = sum / matrix_get(A, i, i);
            }
            const double error = vector_max_diff(x, previous_x);
            if (error < PRECISION) {
                break;
            }
        }
    }
    vector_dealloc(&previous_x);
    return x;
}

// Remember to free the returned vector after calling this function!
Vector gauss_seidel(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = vector_init(b.len, 0.0);
    Vector previous_x = vector_alloc(b.len);
    if (!vector_is_valid(x) || !vector_is_valid(previous_x)) {
        vector_dealloc(&x);
    } else {
        for (size_t k = 0; k < MAX_ITERATIONS; k++) {
            vector_copy_over(&previous_x, x);
            for (size_t i = 0; i < b.len; i++) {
                double sum = b.data[i];
                for (size_t j = 0; j < x.len; j++) {
                    if (i != j) {
                        sum -= matrix_get(A, i, j) * x.data[j];
                    }
                }
                x.data[i] = sum / matrix_get(A, i, i);
            }
            const double error = vector_max_diff(x, previous_x);
            if (error < PRECISION) {
                break;
            }
        }
    }
    vector_dealloc(&previous_x);
    return x;
}

bool columns_condition(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < A.rows; j++) {
            if (i != j) {
                sum += fabs(matrix_get(A, j, i));
            }
        }
        if (fabs(matrix_get(A, i, i)) < sum) {
            return false;
        }
    }
    return true;
}

bool rows_condition(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        double sum = 0.0;
        for (size_t j = 0; j < A.rows; j++) {
            if (i != j) {
                sum += fabs(matrix_get(A, i, j));
            }
        }
        if (fabs(matrix_get(A, i, i)) < sum) {
            return false;
        }
    }
    return true;
}

bool sassenfeld_condition(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    bool condition = false;
    Vector beta = vector_alloc(A.rows);
    if (vector_is_valid(beta)) {
        for (size_t i = 0; i < A.rows; i++) {
            double sum = 0.0;
            for (size_t j = 0; j < i; j++) {
                sum += fabs(matrix_get(A, i, j)) * beta.data[j];
            }
            for (size_t j = (i + 1); j < A.rows; j++) {
                sum += fabs(matrix_get(A, i, j));
            }
            beta.data[i] = sum / matrix_get(A, i, i);
        }
        if (vector_max(beta) < 1.0) {
            condition = true;
        }
        vector_dealloc(&beta);
    }
    return condition;
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