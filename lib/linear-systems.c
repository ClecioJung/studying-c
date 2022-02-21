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
    Vector x = vector_copy(b);
    if (vector_is_valid(x)) {
        back_substitution_over(A, x);
    }
    return x;
}

// Remember to free the returned vector after calling this function!
Vector forward_substitution(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = vector_copy(b);
    if (vector_is_valid(x)) {
        forward_substitution_over(A, x);
    }
    return x;
}

// Remember to free the vector x after calling this function!
System_Type gaussian_elimination(const Matrix A, const Vector b, Vector *const x) {
    if (x == NULL) {
        return Error;
    }
    System_Type type = Error;
    Matrix A_copy = matrix_copy(A);
    *x = vector_copy(b);
    if (matrix_is_valid(A_copy) && vector_is_valid(*x)) {
        type = gaussian_elimination_over(A_copy, *x);
    }
    if (type != Solvable) {
        vector_dealloc(x);
    }
    matrix_dealloc(&A_copy);
    return type;
}

// Remember to free the returned vector after calling this function!
System_Type gauss_jordan(const Matrix A, const Vector b, Vector *const x) {
    if (x == NULL) {
        return Error;
    }
    System_Type type = Error;
    Matrix A_copy = matrix_copy(A);
    *x = vector_copy(b);
    if (matrix_is_valid(A_copy) && vector_is_valid(*x)) {
        type = gauss_jordan_over(A_copy, *x);
    }
    if (type != Solvable) {
        vector_dealloc(x);
    }
    matrix_dealloc(&A_copy);
    return type;
}

// Remember to free the returned vector after calling this function!
System_Type lu_solving(const Matrix A, const Vector b, Vector *const x) {
    if (x == NULL) {
        return Error;
    }
    System_Type type = Error;
    Matrix A_copy = matrix_copy(A);
    *x = vector_copy(b);
    if (matrix_is_valid(A_copy) && vector_is_valid(*x)) {
        type = lu_solving_over(A_copy, *x);
    }
    if (type != Solvable) {
        vector_dealloc(x);
    }
    matrix_dealloc(&A_copy);
    return type;
}

// Remember to free the returned vector after calling this function!
System_Type jacobi_method(const Matrix A, const Vector b, Vector *const x) {
    if (x == NULL) {
        return Error;
    }
    if (A.rows != b.len) {
        return Invalid;
    }
    if (A.cols > A.rows) {
        return Underdetermined;
    }
    if (A.rows > A.cols) {
        // The system may be underdetermined also, in case the extra equations are linearly dependent
        // But this corner case was not implemented here
        return Overdetermined;
    }
    System_Type type = Error;
    *x = vector_init(b.len, 0.0);
    Vector previous_x = vector_alloc(b.len);
    if (vector_is_valid(*x) && vector_is_valid(previous_x)) {
        type = Underdetermined;
        for (size_t k = 0; k < MAX_ITERATIONS; k++) {
            vector_copy_over(previous_x, *x);
            for (size_t i = 0; i < x->len; i++) {
                double sum = b.data[i];
                for (size_t j = 0; j < x->len; j++) {
                    if (i != j) {
                        sum -= matrix_get(A, i, j) * previous_x.data[j];
                    }
                }
                x->data[i] = sum / matrix_get(A, i, i);
            }
            const double error = vector_max_diff(*x, previous_x);
            if (error < PRECISION) {
                type = Solvable;
                break;
            }
        }
    }
    if (type != Solvable) {
        vector_dealloc(x);
    }
    vector_dealloc(&previous_x);
    return type;
}

// Remember to free the returned vector after calling this function!
System_Type gauss_seidel(const Matrix A, const Vector b, Vector *const x) {
    if (x == NULL) {
        return Error;
    }
    if (A.rows != b.len) {
        return Invalid;
    }
    if (A.cols > A.rows) {
        return Underdetermined;
    }
    if (A.rows > A.cols) {
        // The system may be underdetermined also, in case the extra equations are linearly dependent
        // But this corner case was not implemented here
        return Overdetermined;
    }
    System_Type type = Error;
    *x = vector_init(b.len, 0.0);
    Vector previous_x = vector_alloc(b.len);
    if (vector_is_valid(*x) && vector_is_valid(previous_x)) {
        type = Underdetermined;
        for (size_t k = 0; k < MAX_ITERATIONS; k++) {
            vector_copy_over(previous_x, *x);
            for (size_t i = 0; i < b.len; i++) {
                double sum = b.data[i];
                for (size_t j = 0; j < x->len; j++) {
                    if (i != j) {
                        sum -= matrix_get(A, i, j) * x->data[j];
                    }
                }
                x->data[i] = sum / matrix_get(A, i, i);
            }
            if (vector_max_diff(*x, previous_x) < PRECISION) {
                type = Solvable;
                break;
            }
        }
    }
    if (type != Solvable) {
        vector_dealloc(x);
    }
    vector_dealloc(&previous_x);
    return type;
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

// WARNING! This function changes the contents of A and b!
void partial_pivoting_over(const Matrix A, const Vector b) {
    for (size_t k = 0; (k + 1) < A.rows; k++) {
        double w = fabs(matrix_get(A, k, k));
        size_t r = k;
        for (size_t i = (k + 1); i < A.rows; i++) {
            if (fabs(matrix_get(A, i, k)) > w) {
                w = fabs(matrix_get(A, i, k));
                r = i;
            }
        }
        if (r != k) {
            for (size_t i = k; i < A.rows; i++) {
                const double temp = matrix_get(A, k, i);
                matrix_set(A, k, i, matrix_get(A, r, i));
                matrix_set(A, r, i, temp);
            }
            swap(&b.data[k], &b.data[r]);
        }
    }
}

// The vector solution is stored in b
// The matrix is considered to be upper triangular
void back_substitution_over(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return;  // Invalid operation
    }
    for (size_t i = (b.len - 1); i < b.len; i--) {
        for (size_t j = (i + 1); j < b.len; j++) {
            b.data[i] -= matrix_get(A, i, j) * b.data[j];
        }
        b.data[i] /= matrix_get(A, i, i);
    }
}

// The vector solution is stored in b
// The matrix is considered to be lower triangular with ones in the diagonal
void forward_substitution_over(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return;  // Invalid operation
    }
    for (size_t i = 0; i < b.len; i++) {
        for (size_t j = 0; j < i; j++) {
            b.data[i] -= matrix_get(A, i, j) * b.data[j];
        }
    }
}

// WARNING! This function changes the contents of A and b!
// The vector solution is stored in b
System_Type gaussian_elimination_over(const Matrix A, const Vector b) {
    if (A.rows != b.len) {
        return Invalid;
    }
    if (A.cols > A.rows) {
        return Underdetermined;
    }
    if (A.rows > A.cols) {
        // The system may be underdetermined also, in case the extra equations are linearly dependent
        // But this corner case was not implemented here
        return Overdetermined;
    }
    // TODO: check if is a singular system
    partial_pivoting_over(A, b);
    // Forward elimination
    for (size_t k = 0; (k + 1) < A.rows; k++) {
        for (size_t i = (k + 1); i < A.rows; i++) {
            if (are_close(matrix_get(A, k, k), 0.0, PRECISION)) {
                return Underdetermined;
            }
            const double m = matrix_get(A, i, k) / matrix_get(A, k, k);
            for (size_t j = (k + 1); j < A.rows; j++) {
                matrix_dec(A, i, j, m * matrix_get(A, k, j));
            }
            b.data[i] -= m * b.data[k];
        }
    }
    if (are_close(matrix_get(A, (A.rows - 1), (A.rows - 1)), 0.0, PRECISION)) {
        return Underdetermined;
    }
    back_substitution_over(A, b);
    return Solvable;
}

// WARNING! This function changes the contents of A and b!
// The vector solution is stored in b
System_Type gauss_jordan_over(const Matrix A, const Vector b) {
    if (A.rows != b.len) {
        return Invalid;
    }
    if (A.cols > A.rows) {
        return Underdetermined;
    }
    if (A.rows > A.cols) {
        // The system may be underdetermined also, in case the extra equations are linearly dependent
        // But this corner case was not implemented here
        return Overdetermined;
    }
    partial_pivoting_over(A, b);
    for (size_t k = 0; k < A.rows; k++) {
        for (size_t i = 0; i < A.rows; i++) {
            const double p = matrix_get(A, k, k);
            if (are_close(p, 0.0, PRECISION)) {
                return Underdetermined;
            }
            if (i == k) {
                for (size_t j = k; j < A.cols; j++) {
                    matrix_set(A, i, j, matrix_get(A, i, j) / p);
                }
                b.data[i] /= p;
            } else {
                const double m = matrix_get(A, i, k) / p;
                for (size_t j = k; j < A.cols; j++) {
                    matrix_dec(A, i, j, m * matrix_get(A, k, j));
                }
                b.data[i] -= m * b.data[k];
            }
        }
    }
    return Solvable;
}

// WARNING! This function changes the contents of A and b!
// The vector solution is stored in b
System_Type lu_solving_over(const Matrix A, const Vector b) {
    if (A.rows != b.len) {
        return Invalid;
    }
    if (A.cols > A.rows) {
        return Underdetermined;
    }
    if (A.rows > A.cols) {
        // The system may be underdetermined also, in case the extra equations are linearly dependent
        // But this corner case was not implemented here
        return Overdetermined;
    }
    partial_pivoting_over(A, b);
    matrix_lu_dec_over(A);
    // Check if the system is underdetermined
    for (size_t k = 0; k < A.rows; k++) {
        if (are_close(matrix_get(A, k, k), 0.0, PRECISION)) {
            return Underdetermined;
        }
    }
    forward_substitution_over(A, b);
    back_substitution_over(A, b);
    return Solvable;
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