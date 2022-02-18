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

#include "matrix.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "linear-systems.h"
#include "scalar.h"

#define MAX_ITERATIONS 10000
#define PRECISION 1e-9

// Remember to free the matrix after calling this function!
Matrix matrix_alloc(const size_t rows, const size_t cols) {
    Matrix matrix = (Matrix){0};
    if ((rows != 0) && (cols != 0)) {
        matrix.data = (double *)malloc(rows * cols * sizeof(double));
        if (matrix.data != NULL) {
            matrix.rows = rows;
            matrix.cols = cols;
        }
    }
    return matrix;
}

void matrix_dealloc(Matrix *const matrix) {
    if (matrix == NULL) {
        return;
    }
    free(matrix->data);
    matrix->data = NULL;
    matrix->rows = 0;
    matrix->cols = 0;
}

bool matrix_is_valid(const Matrix A) {
    return (A.data != NULL);
}

void matrix_set(const Matrix A, const size_t row, const size_t col, const double value) {
    if (!matrix_is_valid(A) || (row >= A.rows) || (col >= A.cols)) {
        return;
    }
    A.data[row * A.cols + col] = value;
}

double matrix_get(const Matrix A, const size_t row, const size_t col) {
    if (!matrix_is_valid(A) || (row >= A.rows) || (col >= A.cols)) {
        return NAN;
    }
    return A.data[row * A.cols + col];
}

void matrix_inc(const Matrix A, const size_t row, const size_t col, const double value) {
    if (!matrix_is_valid(A) || (row >= A.rows) || (col >= A.cols)) {
        return;
    }
    A.data[row * A.cols + col] += value;
}

void matrix_dec(const Matrix A, const size_t row, const size_t col, const double value) {
    if (!matrix_is_valid(A) || (row >= A.rows) || (col >= A.cols)) {
        return;
    }
    A.data[row * A.cols + col] -= value;
}

bool matrix_is_squared(const Matrix A) {
    return (matrix_is_valid(A) && (A.rows == A.cols));
}

// Remember to free the matrix after calling this function!
Matrix matrix_random(const size_t rows, const size_t cols, const double min, const double max) {
    Matrix matrix = matrix_alloc(rows, cols);
    if (matrix_is_valid(matrix)) {
        srand(time(NULL));
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                matrix_set(matrix, i, j, random_number(min, max));
            }
        }
    }
    return matrix;
}

// Remember to free the matrix after calling this function!
Matrix matrix_init(const size_t rows, const size_t cols, const double value) {
    Matrix matrix = matrix_alloc(rows, cols);
    if (matrix_is_valid(matrix)) {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                matrix_set(matrix, i, j, value);
            }
        }
    }
    return matrix;
}

// Remember to free the matrix after calling this function!
Matrix matrix_identity(const size_t rows) {
    Matrix matrix = matrix_alloc(rows, rows);
    if (matrix_is_valid(matrix)) {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < rows; j++) {
                matrix_set(matrix, i, j, (i == j ? 1.0 : 0.0));
            }
        }
    }
    return matrix;
}

// Remember to free the matrix after calling this function!
Matrix matrix_copy(const Matrix matrix) {
    Matrix new_matrix = matrix_alloc(matrix.rows, matrix.cols);
    if (matrix_is_valid(new_matrix)) {
        for (size_t i = 0; i < matrix.rows; i++) {
            for (size_t j = 0; j < matrix.cols; j++) {
                matrix_set(new_matrix, i, j, matrix_get(matrix, i, j));
            }
        }
    }
    return new_matrix;
}

void matrix_assign(Matrix *const matrix, const Matrix equals) {
    matrix->data = equals.data;
    matrix->rows = equals.rows;
    matrix->cols = equals.cols;
}

void matrix_replace(Matrix *const matrix, const Matrix equals) {
    matrix_dealloc(matrix);
    matrix_assign(matrix, equals);
}

void matrix_print(const Matrix matrix) {
    printf("%*s", 2, "");
    for (size_t j = 0; j < matrix.cols; j++) {
        printf("%*s[%03zu] ", 5, "", j);
    }
    printf("\n");
    for (size_t i = 0; i < matrix.rows; i++) {
        printf("[%03zu]: ", i);
        for (size_t j = 0; j < matrix.cols; j++) {
            const double value = fabs(matrix_get(matrix, i, j)) > PRECISION ? matrix_get(matrix, i, j) : 0.0;
            printf("%-10g ", value);
        }
        printf("\n");
    }
    printf("\n");
}

// Remember to free the returned vector after calling this function!
Vector vector_from_matrix_column(const Matrix A, const size_t col) {
    if (col >= A.cols) {
        return (Vector){0};  // Invalid operation
    }
    Vector vec = vector_alloc(A.rows);
    if (vector_is_valid(vec)) {
        for (size_t i = 0; i < vec.len; i++) {
            vec.data[i] = matrix_get(A, i, col);
        }
    }
    return vec;
}

// Remember to free the returned vector after calling this function!
Matrix matrix_scale(const double scalar, const Matrix A) {
    Matrix result = matrix_alloc(A.rows, A.cols);
    if (matrix_is_valid(result)) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                matrix_set(result, i, j, scalar * matrix_get(A, i, j));
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_sum(const Matrix A, const Matrix B) {
    if ((A.rows != B.rows) || (A.cols != B.cols)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = matrix_alloc(A.rows, A.cols);
    if (matrix_is_valid(result)) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                matrix_set(result, i, j, matrix_get(A, i, j) + matrix_get(B, i, j));
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_sub(const Matrix A, const Matrix B) {
    if ((A.rows != B.rows) || (A.cols != B.cols)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = matrix_alloc(A.rows, A.cols);
    if (matrix_is_valid(result)) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                matrix_set(result, i, j, matrix_get(A, i, j) - matrix_get(B, i, j));
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_mul(const Matrix A, const Matrix B) {
    if (A.cols != B.rows) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = matrix_init(A.rows, B.cols, 0.0);
    if (matrix_is_valid(result)) {
        for (size_t i = 0; i < result.rows; i++) {
            for (size_t j = 0; j < result.cols; j++) {
                for (size_t k = 0; k < A.cols; k++) {
                    matrix_inc(result, i, j, matrix_get(A, i, k) * matrix_get(B, k, j));
                }
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_mul_three(const Matrix A, const Matrix B, const Matrix C) {
    if ((A.cols != B.rows) || (B.cols != C.rows)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = matrix_mul(A, B);
    matrix_replace(&result, matrix_mul(result, C));
    return result;
}

// Remember to free the returned vector after calling this function!
Vector matrix_mul_vector(const Matrix A, const Vector b) {
    if (A.cols != b.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector result = vector_init(A.rows, 0.0);
    if (vector_is_valid(result)) {
        for (size_t i = 0; i < result.len; i++) {
            for (size_t k = 0; k < A.cols; k++) {
                result.data[i] += matrix_get(A, i, k) * b.data[k];
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_transpose(const Matrix A) {
    Matrix transpose = matrix_alloc(A.cols, A.rows);
    if (matrix_is_valid(transpose)) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                matrix_set(transpose, j, i, matrix_get(A, i, j));
            }
        }
    }
    return transpose;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_symmetric(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix sym = matrix_alloc(A.rows, A.cols);
    if (matrix_is_valid(sym)) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                matrix_set(sym, i, j, (matrix_get(A, i, j) + matrix_get(A, j, i)) / 2.0);
            }
        }
    }
    return sym;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_skew_symmetric(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix skew = matrix_alloc(A.rows, A.cols);
    if (matrix_is_valid(skew)) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                matrix_set(skew, i, j, (matrix_get(A, i, j) - matrix_get(A, j, i)) / 2.0);
            }
        }
    }
    return skew;
}

double trace(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return 0.0;  // Invalid operation
    }
    double trace = 0.0;
    for (size_t i = 0; i < A.rows; i++) {
        trace += matrix_get(A, i, i);
    }
    return trace;
}

double determinant(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return 0.0;  // Invalid operation
    }
    // Gaussian elimination
    Matrix B = matrix_copy(A);
    if (!matrix_is_valid(B)) {
        return 0.0;
    }
    for (size_t k = 0; (k + 1) < B.rows; k++) {
        for (size_t i = (k + 1); i < B.rows; i++) {
            const double m = matrix_get(B, i, k) / matrix_get(B, k, k);
            for (size_t j = (k + 1); j < B.rows; j++) {
                matrix_dec(B, i, j, m * matrix_get(B, k, j));
            }
        }
    }
    double det = 1.0;
    for (size_t k = 0; k < B.rows; k++) {
        det *= matrix_get(B, k, k);
    }
    matrix_dealloc(&B);
    return det;
}

// Remember to free L and U matrices after calling this function!
void lu_decomposition(const Matrix A, Matrix *const L, Matrix *const U) {
    if (!matrix_is_squared(A) || (L == NULL) || (U == NULL)) {
        return;  // Invalid operation
    }
    *L = matrix_identity(A.rows);
    *U = matrix_copy(A);
    if (!matrix_is_valid(*L) || !matrix_is_valid(*U)) {
        matrix_dealloc(L);
        matrix_dealloc(U);
        return;
    }
    for (size_t k = 0; (k + 1) < U->rows; k++) {
        for (size_t i = (k + 1); i < U->rows; i++) {
            matrix_set(*L, i, k, matrix_get(*U, i, k) / matrix_get(*U, k, k));
            matrix_set(*U, i, k, 0.0);
            for (size_t j = (k + 1); j < U->rows; j++) {
                matrix_dec(*U, i, j, matrix_get(*L, i, k) * matrix_get(*U, k, j));
            }
        }
    }
}

// Remember to free L and U matrices after calling this function!
void lu_crout_decomposition(const Matrix A, Matrix *const L, Matrix *const U) {
    if (!matrix_is_squared(A) || (L == NULL) || (U == NULL)) {
        return;  // Invalid operation
    }
    // Crout Decomposition
    *L = matrix_init(A.rows, A.cols, 0.0);
    *U = matrix_identity(A.rows);
    if (!matrix_is_valid(*L) || !matrix_is_valid(*U)) {
        matrix_dealloc(L);
        matrix_dealloc(U);
        return;
    }
    for (size_t k = 0; k < A.rows; k++) {
        for (size_t i = k; i < A.rows; i++) {
            matrix_set(*L, i, k, matrix_get(A, i, k));
            for (size_t l = 0; l < k; l++) {
                matrix_dec(*L, i, k, matrix_get(*L, i, l) * matrix_get(*U, l, k));
            }
        }
        if ((k + 1) != A.rows) {
            for (size_t j = (k + 1); j < A.rows; j++) {
                matrix_set(*U, k, j, matrix_get(A, k, j));
                for (size_t l = 0; l < k; l++) {
                    matrix_dec(*U, k, j, matrix_get(*L, k, l) * matrix_get(*U, l, j));
                }
                matrix_set(*U, k, j, matrix_get(*U, k, j) / matrix_get(*L, k, k));
            }
        }
    }
}

// Remember to free Q and R matrices after calling this function!
void qr_decomposition(const Matrix A, Matrix *const Q, Matrix *const R) {
    if (!matrix_is_squared(A) || (Q == NULL) || (R == NULL)) {
        return;  // Invalid operation
    }
    *Q = matrix_copy(A);
    *R = matrix_init(A.rows, A.cols, 0.0);
    Vector a = vector_alloc(A.rows);
    Vector q = vector_alloc(A.rows);
    if (!matrix_is_valid(*Q) || !matrix_is_valid(*R) || !vector_is_valid(a) || !vector_is_valid(q)) {
        matrix_dealloc(Q);
        matrix_dealloc(R);
    } else {
        for (size_t i = 0; i < A.cols; i++) {
            vector_from_matrix_column_over(&a, A, i);
            // Gram-Schmidt process
            for (size_t j = 0; j < i; j++) {
                vector_from_matrix_column_over(&q, *Q, j);
                matrix_set(*R, j, i, dot_product(q, a));
                for (size_t k = 0; k < A.rows; k++) {
                    matrix_dec(*Q, k, i, matrix_get(*R, j, i) * q.data[k]);
                }
            }
            vector_from_matrix_column_over(&q, *Q, i);
            matrix_set(*R, i, i, euclidean_norm(q));
            for (size_t k = 0; k < A.rows; k++) {
                matrix_set(*Q, k, i, q.data[k] / matrix_get(*R, i, i));
            }
        }
    }
    vector_dealloc(&a);
    vector_dealloc(&q);
}

// Remember to free the returned matrix after calling this function!
Matrix householder_matrix(const Vector vec) {
    Matrix matrix = matrix_alloc(vec.len, vec.len);
    if (matrix_is_valid(matrix)) {
        double norm_squared = 0.0;
        for (size_t i = 0; i < vec.len; i++) {
            norm_squared += vec.data[i] * vec.data[i];
        }
        for (size_t i = 0; i < matrix.rows; i++) {
            for (size_t j = 0; j < matrix.cols; j++) {
                matrix_set(matrix, i, j, ((i == j) ? 1.0 : 0.0) - (2.0 * vec.data[i] * vec.data[j]) / norm_squared);
            }
        }
    }
    return matrix;
}

// Remember to free the matrices U and H after calling this function!
// Decomposition A = U * H * U^T,
// with H in upper Hessenberg form and U orthogonal
void upper_hessenberg_matrix(const Matrix A, Matrix *const U, Matrix *const H) {
    if (!matrix_is_squared(A) || (U == NULL) || (H == NULL)) {
        return;  // Invalid operation
    }
    *U = matrix_identity(A.rows);
    *H = matrix_copy(A);
    Matrix P = matrix_alloc(A.rows, A.cols);
    Vector v = vector_alloc(A.rows);
    if (!matrix_is_valid(*H) || !matrix_is_valid(*U) || !matrix_is_valid(P) || !vector_is_valid(v)) {
        matrix_dealloc(U);
        matrix_dealloc(H);
    } else {
        for (size_t i = 1; (i + 1) < A.rows; i++) {
            // Calculate Householder matrix
            for (size_t j = 0; j < v.len; j++) {
                v.data[j] = (j >= i) ? matrix_get(*H, j, i - 1) : 0.0;
            }
            v.data[i] += sign(v.data[i]) * euclidean_norm(v);
            double norm_squared = 0.0;
            for (size_t i = 0; i < v.len; i++) {
                norm_squared += v.data[i] * v.data[i];
            }
            for (size_t i = 0; i < P.rows; i++) {
                for (size_t j = 0; j < P.cols; j++) {
                    matrix_set(P, i, j, ((i == j) ? 1.0 : 0.0) - (2.0 * v.data[i] * v.data[j]) / norm_squared);
                }
            }
            matrix_replace(H, matrix_mul_three(P, *H, P));
            matrix_replace(U, matrix_mul(*U, P));
            if (!matrix_is_valid(*H) || !matrix_is_valid(*U)) {
                matrix_dealloc(U);
                matrix_dealloc(H);
                break;
            }
        }
    }
    matrix_dealloc(&P);
    vector_dealloc(&v);
}

// Remember to free the matrices U and T after calling this function!
void schur_decomposition(const Matrix A, Matrix *const U, Matrix *const T) {
    if (!matrix_is_squared(A) || (U == NULL) || (T == NULL)) {
        return;  // Invalid operation
    }
    upper_hessenberg_matrix(A, U, T);
    if (!matrix_is_valid(*U) || !matrix_is_valid(*T)) {
    schur_decomposition_safe_exit:
        matrix_dealloc(U);
        matrix_dealloc(T);
        return;
    }
    // QR algorithm
    for (size_t i = 0; i < MAX_ITERATIONS; i++) {
        Matrix Q = (Matrix){0};
        Matrix R = (Matrix){0};
        qr_decomposition(*T, &Q, &R);
        if (!matrix_is_valid(R) || !matrix_is_valid(Q)) {
            matrix_dealloc(&Q);
            matrix_dealloc(&R);
            goto schur_decomposition_safe_exit;
        }
        matrix_replace(T, matrix_mul(R, Q));
        matrix_replace(U, matrix_mul(*U, Q));
        matrix_dealloc(&Q);
        matrix_dealloc(&R);
        if (!matrix_is_valid(*U) || !matrix_is_valid(*T)) {
            goto schur_decomposition_safe_exit;
        }
        if (matrix_is_upper_triangular(*T)) {
            break;
        }
    }
}

// Remember to free the returned vector after calling this function!
Vector eigenvalues(const Matrix A) {
    Matrix U = (Matrix){0};
    Matrix T = (Matrix){0};
    schur_decomposition(A, &U, &T);
    Vector eig = vector_alloc(T.rows);
    if (vector_is_valid(eig) && matrix_is_upper_triangular(T)) {
        for (size_t i = 0; i < T.rows; i++) {
            eig.data[i] = matrix_get(T, i, i);
        }
    }
    matrix_dealloc(&U);
    matrix_dealloc(&T);
    return eig;
}

// Remember to free vec after calling this function!
double power_method(const Matrix A, Vector *const vec) {
    if (!matrix_is_squared(A) || (vec == NULL)) {
        return NAN;  // Invalid operation
    }
    *vec = vector_random(A.rows, 0.0, 1.0);
    Vector previous_vec = vector_alloc(vec->len);
    if (!vector_is_valid(*vec) || !vector_is_valid(previous_vec)) {
    power_method_safe_exit:
        vector_dealloc(vec);
        vector_dealloc(&previous_vec);
        return NAN;
    }
    double eig = 0.0;
    for (size_t i = 0; i < MAX_ITERATIONS; i++) {
        vector_copy_over(&previous_vec, *vec);
        vector_replace(vec, matrix_mul_vector(A, *vec));
        if (!vector_is_valid(*vec)) {
            goto power_method_safe_exit;
        }
        eig = euclidean_norm(*vec);
        vector_scale_over((1.0 / eig), vec);
        if (vector_max_diff(*vec, previous_vec) < PRECISION) {
            break;
        }
    }
    vector_dealloc(&previous_vec);
    return eig;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_inverse(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return (Matrix){0};  // Invalid operation
    }
    // Computes the inverse matrix of A using LU decomposition
    Matrix L, U;
    lu_decomposition(A, &L, &U);
    Matrix inv = matrix_alloc(A.rows, A.cols);
    Vector b = vector_alloc(A.rows);
    if (!matrix_is_valid(L) || !matrix_is_valid(U) || !matrix_is_valid(inv) || !vector_is_valid(b)) {
        matrix_dealloc(&L);
        matrix_dealloc(&U);
        matrix_dealloc(&inv);
        vector_dealloc(&b);
        return (Matrix){0};
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.rows; j++) {
            b.data[j] = (j == i) ? 1.0 : 0.0;
        }
        Vector d = forward_substitution(L, b);
        if (!vector_is_valid(d)) {
            matrix_dealloc(&inv);
            break;
        }
        Vector x = back_substitution(U, d);
        if (!vector_is_valid(x)) {
            vector_dealloc(&d);
            matrix_dealloc(&inv);
            break;
        }
        for (size_t j = 0; j < A.rows; j++) {
            matrix_set(inv, j, i, x.data[j]);
        }
        vector_dealloc(&d);
        vector_dealloc(&x);
    }
    matrix_dealloc(&L);
    matrix_dealloc(&U);
    vector_dealloc(&b);
    return inv;
}

Matrix pseudo_inverse(const Matrix A) {
    Matrix temp;
    Matrix transpose = matrix_transpose(A);
    if (A.rows > A.cols) {
        // Left pseudo inverse
        // A^+ = (A^T * A)^-1 * A^T
        temp = matrix_mul(transpose, A);
        matrix_replace(&temp, matrix_inverse(temp));
        matrix_replace(&temp, matrix_mul(temp, transpose));
    } else {
        // Right pseudo inverse
        // A^+ = A^T * (A * A^T)^-1
        temp = matrix_mul(A, transpose);
        matrix_replace(&temp, matrix_inverse(temp));
        matrix_replace(&temp, matrix_mul(transpose, temp));
    }
    matrix_dealloc(&transpose);
    return temp;
}

double matrix_max_diff(const Matrix A, const Matrix previous_A) {
    if ((A.rows != previous_A.rows) || (A.cols != previous_A.cols)) {
        return INFINITY;  // Invalid operation
    }
    double error = 0.0;
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            error = maximum(fabs(matrix_get(A, i, j) - matrix_get(previous_A, i, j)), error);
        }
    }
    return error;
}

bool matrix_are_equal(const Matrix A, const Matrix B) {
    if ((A.rows != B.rows) || (A.cols != B.cols)) {
        return false;
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            if (!are_close(matrix_get(A, i, j), matrix_get(B, i, j), PRECISION)) {
                return false;
            }
        }
    }
    return true;
}

bool matrix_is_orthogonal(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;
    }
    Matrix temp = matrix_transpose(A);
    matrix_replace(&temp, matrix_mul(temp, A));
    Matrix identity = matrix_identity(A.rows);
    const bool are_equal = matrix_are_equal(temp, identity);
    matrix_dealloc(&temp);
    matrix_dealloc(&identity);
    return are_equal;
}

bool matrix_is_upper_triangular(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    double error = 0.0;
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < i; j++) {
            error = maximum(error, fabs(matrix_get(A, i, j)));
        }
    }
    return (error < PRECISION);
}

void vector_from_matrix_column_over(const Vector *const vector, const Matrix A, const size_t col) {
    if ((vector->len != A.rows) || (col >= A.cols)) {
        return;  // Invalid operation
    }
    for (size_t i = 0; i < vector->len; i++) {
        vector->data[i] = matrix_get(A, i, col);
    }
}

void matrix_identity_over(const Matrix *const matrix) {
    if (matrix_is_squared(*matrix)) {
        for (size_t i = 0; i < matrix->rows; i++) {
            for (size_t j = 0; j < matrix->rows; j++) {
                matrix_set(*matrix, i, j, (i == j ? 1.0 : 0.0));
            }
        }
    }
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