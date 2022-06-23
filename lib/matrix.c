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
#define ITERATION_PRECISION 1e-10
#define COMPARATION_PRECISION 1e-9

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
        matrix_init_over(matrix, value);
    }
    return matrix;
}

// Remember to free the matrix after calling this function!
Matrix matrix_identity(const size_t rows) {
    Matrix matrix = matrix_alloc(rows, rows);
    if (matrix_is_valid(matrix)) {
        matrix_identity_over(matrix);
    }
    return matrix;
}

// Remember to free the matrix after calling this function!
Matrix matrix_copy(const Matrix matrix) {
    Matrix new_matrix = matrix_alloc(matrix.rows, matrix.cols);
    if (matrix_is_valid(new_matrix)) {
        matrix_copy_over(new_matrix, matrix);
    }
    return new_matrix;
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
            const double value = fabs(matrix_get(matrix, i, j)) > COMPARATION_PRECISION ? matrix_get(matrix, i, j) : 0.0;
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
        vector_from_matrix_column_over(vec, A, col);
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
        matrix_sum_over(result, A, B);
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
        matrix_sub_over(result, A, B);
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_mul(const Matrix A, const Matrix B) {
    if (A.cols != B.rows) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = matrix_alloc(A.rows, B.cols);
    if (matrix_is_valid(result)) {
        matrix_mul_over(result, A, B);
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_mul_three(const Matrix A, const Matrix B, const Matrix C) {
    if ((A.cols != B.rows) || (B.cols != C.rows)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix mul = matrix_mul(A, B);
    Matrix result = matrix_mul(mul, C);
    matrix_dealloc(&mul);
    return result;
}

// Remember to free the returned vector after calling this function!
Vector matrix_mul_vector(const Matrix A, const Vector b) {
    if (A.cols != b.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector result = vector_alloc(A.rows);
    if (vector_is_valid(result)) {
        matrix_mul_vector_over(result, A, b);
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

double matrix_trace(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return NAN;  // Invalid operation
    }
    double trace = 0.0;
    for (size_t i = 0; i < A.rows; i++) {
        trace += matrix_get(A, i, i);
    }
    return trace;
}

double matrix_determinant(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return NAN;  // Invalid operation
    }
    Matrix A_copy = matrix_copy(A);
    if (!matrix_is_valid(A_copy)) {
        return NAN;
    }
    // Forward elimination (Gauss)
    for (size_t k = 0; (k + 1) < A_copy.rows; k++) {
        for (size_t i = (k + 1); i < A_copy.rows; i++) {
            const double m = matrix_get(A_copy, i, k) / matrix_get(A_copy, k, k);
            for (size_t j = (k + 1); j < A_copy.rows; j++) {
                matrix_dec(A_copy, i, j, m * matrix_get(A_copy, k, j));
            }
        }
    }
    double det = 1.0;
    for (size_t k = 0; k < A_copy.rows; k++) {
        det *= matrix_get(A_copy, k, k);
    }
    matrix_dealloc(&A_copy);
    return det;
}

// Remember to free L and U matrices after calling this function!
void matrix_lu_decomposition(const Matrix A, Matrix *const L, Matrix *const U) {
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
void matrix_lu_crout_decomposition(const Matrix A, Matrix *const L, Matrix *const U) {
    if (!matrix_is_squared(A) || (L == NULL) || (U == NULL)) {
        return;  // Invalid operation
    }
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
void matrix_qr_decomposition(const Matrix A, Matrix *const Q, Matrix *const R) {
    if (!matrix_is_squared(A) || (Q == NULL) || (R == NULL)) {
        return;  // Invalid operation
    }
    *Q = matrix_copy(A);
    *R = matrix_alloc(A.rows, A.cols);
    if (!matrix_is_valid(*Q) || !matrix_is_valid(*R)) {
        matrix_dealloc(Q);
        matrix_dealloc(R);
    } else {
        matrix_qr_dec_over(*Q, *R);
    }
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_householder(const Vector vec) {
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
void matrix_upper_hessenberg(const Matrix A, Matrix *const H, Matrix *const U) {
    if (!matrix_is_squared(A) || (U == NULL) || (H == NULL)) {
        return;  // Invalid operation
    }
    *U = matrix_alloc(A.rows, A.cols);
    *H = matrix_copy(A);
    if (!matrix_is_valid(*H) || !matrix_is_valid(*U)) {
        matrix_dealloc(U);
        matrix_dealloc(H);
    } else {
        matrix_upper_hessenberg_over(*H, *U);
    }
}

// Remember to free the matrices U and T after calling this function!
// Decomposition A = U * T * U^T,
// with T in upper triangular form and U orthogonal
void matrix_schur_decomposition(const Matrix A, Matrix *const T, Matrix *const U) {
    if (!matrix_is_squared(A) || (U == NULL) || (T == NULL)) {
        return;  // Invalid operation
    }
    *T = matrix_copy(A);
    *U = matrix_alloc(A.rows, A.cols);
    if (!matrix_is_valid(*T) || !matrix_is_valid(*U)) {
        matrix_dealloc(T);
        matrix_dealloc(U);
    } else {
        matrix_schur_dec_over(*T, *U);
    }
}

// Matrix diagonalization
// Remember to free P and D matrices after calling this function!
// Decomposition A = P * D * P^-1, with D in diagonal form
// Only work for real non repeated eigenvalues
// Based on: https://math.stackexchange.com/questions/3947108/how-to-get-eigenvectors-using-qr-algorithm
void matrix_diagonalization(const Matrix A, Matrix *const P, Matrix *const D) {
    if (!matrix_is_squared(A) || (P == NULL) || (D == NULL)) {
        return;  // Invalid operation
    }
    *D = matrix_copy(A);
    *P = matrix_alloc(A.rows, A.cols);
    if (!matrix_is_valid(*D) || !matrix_is_valid(*P)) {
        matrix_dealloc(D);
        matrix_dealloc(P);
    } else {
        matrix_diagonalization_over(*D, *P);
    }
}

// Remember to free the returned vector after calling this function!
Vector matrix_eigenvalues(const Matrix A) {
    Matrix U = (Matrix){0};
    Matrix T = (Matrix){0};
    matrix_schur_decomposition(A, &T, &U);
    Vector eig = vector_alloc(T.rows);
    if (vector_is_valid(eig) && matrix_is_upper_triangular(T)) {
        for (size_t i = 0; i < T.rows; i++) {
            eig.data[i] = matrix_get(T, i, i);
        }
    } else {
        vector_dealloc(&eig);
    }
    matrix_dealloc(&U);
    matrix_dealloc(&T);
    return eig;
}

// Remember to free vec after calling this function!
double matrix_power_method(const Matrix A, Vector *const vec) {
    if (!matrix_is_squared(A) || (vec == NULL)) {
        return NAN;  // Invalid operation
    }
    *vec = vector_random(A.rows, 0.0, 1.0);
    Vector previous_vec = vector_alloc(vec->len);
    if (!vector_is_valid(*vec) || !vector_is_valid(previous_vec)) {
        vector_dealloc(vec);
        vector_dealloc(&previous_vec);
        return NAN;
    }
    double eig = 0.0;
    for (size_t i = 0; i < MAX_ITERATIONS; i++) {
        vector_copy_over(previous_vec, *vec);
        matrix_mul_vector_over(*vec, A, previous_vec);
        eig = vector_norm(*vec);
        vector_scale_over((1.0 / eig), *vec);
        if (vector_max_diff(*vec, previous_vec) < ITERATION_PRECISION) {
            break;
        }
    }
    vector_dealloc(&previous_vec);
    return eig;
}

// Remember to free vec after calling this function!
double matrix_inverse_power_method(const Matrix A, Vector *const vec) {
    if (!matrix_is_squared(A) || (vec == NULL)) {
        return NAN;  // Invalid operation
    }
    Matrix invA = matrix_inverse(A);
    *vec = vector_random(A.rows, 0.0, 1.0);
    Vector previous_vec = vector_alloc(vec->len);
    if (!vector_is_valid(*vec) || !vector_is_valid(previous_vec)) {
        vector_dealloc(vec);
        vector_dealloc(&previous_vec);
        return NAN;
    }
    double eig = 0.0;
    for (size_t i = 0; i < MAX_ITERATIONS; i++) {
        vector_copy_over(previous_vec, *vec);
        matrix_mul_vector_over(*vec, invA, previous_vec);
        eig = -1.0 / vector_norm(*vec);
        vector_scale_over(eig, *vec);
        if (vector_max_diff(*vec, previous_vec) < ITERATION_PRECISION) {
            break;
        }
    }
    matrix_dealloc(&invA);
    vector_dealloc(&previous_vec);
    return eig;
}

// Remember to free the returned matrix after calling this function!
// Decomposition: A = L * L^T, with L lower-triangular.
// This decomposition is only possible for symmetric matrices
Matrix matrix_cholesky_decomposition(const Matrix A) {
    if (!matrix_is_symmetric(A)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix L = matrix_init(A.rows, A.cols, 0.0);
    if (matrix_is_valid(L)) {
        for (size_t i = 0; i < L.rows; i++) {
            for (size_t j = 0; j <= i; j++) {
                double sum = 0.0;
                for (size_t k = 0; k < j; k++) {
                    sum += matrix_get(L, i, k) * matrix_get(L, j, k);
                }
                if (i == j) {
                    matrix_set(L, i, i, square_root(matrix_get(A, i, i) - sum));
                } else {
                    matrix_set(L, i, j, (matrix_get(A, i, j) - sum) / matrix_get(L, j, j));
                }
            }
        }
    }
    return L;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_inverse(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix A_copy = matrix_copy(A);
    Matrix inv = matrix_alloc(A.rows, A.cols);
    Vector b = vector_alloc(A.rows);
    if (!matrix_is_valid(A_copy) || !matrix_is_valid(inv) || !vector_is_valid(b)) {
        matrix_dealloc(&inv);
    } else {
        // Computes the inverse matrix of A using LU decomposition
        matrix_lu_dec_over(A_copy);
        for (size_t i = 0; i < A_copy.rows; i++) {
            for (size_t j = 0; j < A_copy.rows; j++) {
                b.data[j] = (j == i) ? 1.0 : 0.0;
            }
            forward_substitution_over(A_copy, b);
            back_substitution_over(A_copy, b);
            for (size_t j = 0; j < A.rows; j++) {
                matrix_set(inv, j, i, b.data[j]);
            }
        }
    }
    matrix_dealloc(&A_copy);
    vector_dealloc(&b);
    return inv;
}

Matrix matrix_pseudo_inverse(const Matrix A) {
    Matrix pseudo_inv, mul, inv;
    Matrix transpose = matrix_transpose(A);
    if (A.rows > A.cols) {
        // Left pseudo inverse
        // A^+ = (A^T * A)^-1 * A^T
        mul = matrix_mul(transpose, A);
        inv = matrix_inverse(mul);
        pseudo_inv = matrix_mul(inv, transpose);
    } else {
        // Right pseudo inverse
        // A^+ = A^T * (A * A^T)^-1
        mul = matrix_mul(A, transpose);
        inv = matrix_inverse(mul);
        pseudo_inv = matrix_mul(transpose, inv);
    }
    matrix_dealloc(&mul);
    matrix_dealloc(&inv);
    matrix_dealloc(&transpose);
    return pseudo_inv;
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
            if (!are_close(matrix_get(A, i, j), matrix_get(B, i, j), COMPARATION_PRECISION)) {
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
    Matrix transpose = matrix_transpose(A);
    Matrix mul = matrix_mul(transpose, A);
    Matrix identity = matrix_identity(A.rows);
    const bool are_equal = matrix_are_equal(mul, identity);
    matrix_dealloc(&transpose);
    matrix_dealloc(&mul);
    matrix_dealloc(&identity);
    return are_equal;
}

bool matrix_is_diagonal(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            if (i != j) {
                if (fabs(matrix_get(A, i, j)) >= COMPARATION_PRECISION) {
                    return false;
                }
            }
        }
    }
    return true;
}

bool matrix_is_upper_triangular(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < i; j++) {
            if (fabs(matrix_get(A, i, j)) >= COMPARATION_PRECISION) {
                return false;
            }
        }
    }
    return true;
}

bool matrix_is_lower_triangular(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = (i + 1); j < A.cols; j++) {
            if (fabs(matrix_get(A, i, j)) >= COMPARATION_PRECISION) {
                return false;
            }
        }
    }
    return true;
}

bool matrix_is_symmetric(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < i; j++) {
            if (!are_close(matrix_get(A, i, j), matrix_get(A, j, i), COMPARATION_PRECISION)) {
                return false;
            }
        }
    }
    return true;
}

bool matrix_is_skew_symmetric(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j <= i; j++) {
            if (!are_close(matrix_get(A, i, j), -matrix_get(A, j, i), COMPARATION_PRECISION)) {
                return false;
            }
        }
    }
    return true;
}

bool matrix_is_null_space(const Matrix A, const Vector vec) {
    Vector result = matrix_mul_vector(A, vec);
    const bool is_null_space = vector_is_null(result);
    vector_dealloc(&result);
    return is_null_space;
}

void matrix_init_over(const Matrix A, const double value) {
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            matrix_set(A, i, j, value);
        }
    }
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_jacobian(const multivariable_fn fn, const Vector x) {
    Vector y = fn(x);
    Matrix J = matrix_alloc(y.len, x.len);
    if ((!vector_is_valid(y)) || (!matrix_is_valid(J))) {
        vector_dealloc(&y);
        matrix_dealloc(&J);
        return (Matrix){0};
    }
    const double increment = ITERATION_PRECISION;
    // Computes the Jacobian matrix J
    for (size_t j = 0; j < J.cols; j++) {
        // Increment x_j in order to calculate derivative
        x.data[j] += increment;
        Vector yj = fn(x);
        if (!vector_is_valid(yj)) {
            x.data[j] -= increment;  // Undo increment
            vector_dealloc(&y);
            matrix_dealloc(&J);
            return (Matrix){0};
        }
        for (size_t i = 0; i < J.rows; i++) {
            const double diff = (yj.data[i] - y.data[i]) / increment;
            matrix_set(J, i, j, diff);
        }
        x.data[j] -= increment;  // Undo increment
        vector_dealloc(&yj);
    }
    vector_dealloc(&y);
    return J;
}

void vector_from_matrix_column_over(const Vector vector, const Matrix A, const size_t col) {
    if ((vector.len != A.rows) || (col >= A.cols)) {
        return;  // Invalid operation
    }
    for (size_t i = 0; i < vector.len; i++) {
        vector.data[i] = matrix_get(A, i, col);
    }
}

void matrix_scale_over(const double scalar, const Matrix A) {
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            matrix_set(A, i, j, scalar * matrix_get(A, i, j));
        }
    }
}

void matrix_sum_over(const Matrix result, const Matrix A, const Matrix B) {
    if ((result.rows != A.rows) || (result.cols != A.cols) || (A.rows != B.rows) || (A.cols != B.cols)) {
        return;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            const double value = matrix_get(A, i, j) + matrix_get(B, i, j);
            matrix_set(result, i, j, value);
        }
    }
}

void matrix_sub_over(const Matrix result, const Matrix A, const Matrix B) {
    if ((result.rows != A.rows) || (result.cols != A.cols) || (A.rows != B.rows) || (A.cols != B.cols)) {
        return;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            const double value = matrix_get(A, i, j) - matrix_get(B, i, j);
            matrix_set(result, i, j, value);
        }
    }
}

void matrix_identity_over(const Matrix matrix) {
    if (matrix_is_squared(matrix)) {
        for (size_t i = 0; i < matrix.rows; i++) {
            for (size_t j = 0; j < matrix.rows; j++) {
                matrix_set(matrix, i, j, (i == j ? 1.0 : 0.0));
            }
        }
    }
}

void matrix_copy_over(const Matrix matrix, const Matrix to_copy) {
    if ((matrix.rows != to_copy.rows) || (matrix.cols != to_copy.cols)) {
        return;  // Invalid operation
    }
    for (size_t i = 0; i < matrix.rows; i++) {
        for (size_t j = 0; j < matrix.cols; j++) {
            matrix_set(matrix, i, j, matrix_get(to_copy, i, j));
        }
    }
}

// This function stores the result of the multiplication A * B in the Matrix result.
// This matrix must be preallocated before calling this function!
void matrix_mul_over(const Matrix result, const Matrix A, const Matrix B) {
    if ((A.cols != B.rows) || (result.rows != A.rows) || (result.cols != B.cols)) {
        return;  // Invalid operation
    }
    matrix_init_over(result, 0.0);
    for (size_t i = 0; i < result.rows; i++) {
        for (size_t j = 0; j < result.cols; j++) {
            for (size_t k = 0; k < A.cols; k++) {
                matrix_inc(result, i, j, matrix_get(A, i, k) * matrix_get(B, k, j));
            }
        }
    }
}

// This function stores the result of the multiplication A * b in the Vector result.
// This vector must be preallocated before calling this function!
void matrix_mul_vector_over(const Vector result, const Matrix A, const Vector b) {
    if ((A.cols != b.len) || (result.len != A.rows)) {
        return;  // Invalid operation
    }
    vector_init_over(result, 0.0);
    for (size_t i = 0; i < result.len; i++) {
        for (size_t k = 0; k < A.cols; k++) {
            result.data[i] += matrix_get(A, i, k) * b.data[k];
        }
    }
}

// LU decomposition (version with less memory usage)
// In this version, the LU matrices are stored in the matrix A, changing its contents!
void matrix_lu_dec_over(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return;  // Invalid operation
    }
    for (size_t k = 0; (k + 1) < A.rows; k++) {
        for (size_t i = (k + 1); i < A.rows; i++) {
            matrix_set(A, i, k, matrix_get(A, i, k) / matrix_get(A, k, k));
            for (size_t j = (k + 1); j < A.rows; j++) {
                matrix_dec(A, i, j, matrix_get(A, i, k) * matrix_get(A, k, j));
            }
        }
    }
}

// LU Crout decomposition (version with less memory usage)
// In this version, the LU matrices are stored in the matrix A, changing its contents!
void matrix_lu_crout_dec_over(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return;  // Invalid operation
    }
    for (size_t k = 0; k < A.rows; k++) {
        for (size_t i = k; i < A.rows; i++) {
            for (size_t l = 0; l < k; l++) {
                matrix_dec(A, i, k, matrix_get(A, i, l) * matrix_get(A, l, k));
            }
        }
        if ((k + 1) != A.rows) {
            for (size_t j = (k + 1); j < A.rows; j++) {
                for (size_t l = 0; l < k; l++) {
                    matrix_dec(A, k, j, matrix_get(A, k, l) * matrix_get(A, l, j));
                }
                matrix_set(A, k, j, matrix_get(A, k, j) / matrix_get(A, k, k));
            }
        }
    }
}

// This function undo the effect of lu_dec_over()
// This is equivalent of performing the multiplication A = L*U
// So, the contents of the matrix A are sobrescribed
void matrix_undo_lu_over(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return;  // Invalid operation
    }
    for (size_t j = (A.cols - 1); j < A.cols; j--) {
        for (size_t i = (A.rows - 1); i < A.rows; i--) {
            if (i <= j) {
                for (size_t k = 0; k < i; k++) {
                    matrix_inc(A, i, j, matrix_get(A, i, k) * matrix_get(A, k, j));
                }
            } else {
                matrix_set(A, i, j, matrix_get(A, i, j) * matrix_get(A, j, j));
                for (size_t k = 0; k < j; k++) {
                    matrix_inc(A, i, j, matrix_get(A, i, k) * matrix_get(A, k, j));
                }
            }
        }
    }
}

// This function undo the effect of lu_dec_crout_over()
// This is equivalent of performing the multiplication A = L*U
// So, the contents of the matrix A are sobrescribed
void matrix_undo_lu_crout_over(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return;  // Invalid operation
    }
    for (size_t j = (A.cols - 1); j < A.cols; j--) {
        for (size_t i = (A.rows - 1); i < A.rows; i--) {
            if (i < j) {
                matrix_set(A, i, j, matrix_get(A, i, j) * matrix_get(A, i, i));
            }
            if (i <= j) {
                for (size_t k = 0; k < i; k++) {
                    matrix_inc(A, i, j, matrix_get(A, i, k) * matrix_get(A, k, j));
                }
            } else {
                for (size_t k = 0; k < j; k++) {
                    matrix_inc(A, i, j, matrix_get(A, i, k) * matrix_get(A, k, j));
                }
            }
        }
    }
}

// QR decomposition (version with less memory usage)
// In this version, the Q matrix is stored in the matrix A, changing its contents!
// The matrix R must be preallocated!
// WARNING! This function changes the contents of A and R!
void matrix_qr_dec_over(const Matrix A, const Matrix R) {
    if (!matrix_is_squared(A) || !matrix_is_squared(R) || (A.rows != R.rows)) {
        return;  // Invalid operation
    }
    matrix_init_over(R, 0.0);
    for (size_t i = 0; i < A.cols; i++) {
        // Gram-Schmidt process
        for (size_t j = 0; j < i; j++) {
            double dot_product_ai_aj = 0.0;
            for (size_t k = 0; k < A.rows; k++) {
                dot_product_ai_aj += matrix_get(A, k, i) * matrix_get(A, k, j);
            }
            matrix_set(R, j, i, dot_product_ai_aj);
            for (size_t k = 0; k < A.rows; k++) {
                matrix_dec(A, k, i, dot_product_ai_aj * matrix_get(A, k, j));
            }
        }
        double norm_ai = 0.0;
        for (size_t k = 0; k < A.rows; k++) {
            norm_ai += matrix_get(A, k, i) * matrix_get(A, k, i);
        }
        norm_ai = square_root(norm_ai);
        matrix_set(R, i, i, norm_ai);
        for (size_t k = 0; k < A.rows; k++) {
            matrix_set(A, k, i, matrix_get(A, k, i) / norm_ai);
        }
    }
}

// upper Hessenberg matrix decomposition (version with less memory usage)
// In this version, the H matrix is stored in the matrix A, changing its contents!
// The matrix U must be preallocated!
// WARNING! This function changes the contents of A and U!
// Decomposition A = U * H * U^T,
// with H in upper Hessenberg form and U orthogonal
void matrix_upper_hessenberg_over(const Matrix A, const Matrix U) {
    if (!matrix_is_squared(A) || !matrix_is_squared(U) || (A.rows != U.rows)) {
        return;  // Invalid operation
    }
    Matrix P = matrix_alloc(A.rows, A.cols);
    Matrix temp = matrix_alloc(A.rows, A.cols);
    Vector v = vector_alloc(A.rows);
    if (matrix_is_valid(P) && vector_is_valid(v) && matrix_is_valid(temp)) {
        matrix_identity_over(U);
        for (size_t i = 1; (i + 1) < A.rows; i++) {
            // Calculate Householder matrix
            for (size_t j = 0; j < v.len; j++) {
                v.data[j] = (j >= i) ? matrix_get(A, j, i - 1) : 0.0;
            }
            v.data[i] += sign(v.data[i]) * vector_norm(v);
            double norm_squared = 0.0;
            for (size_t i = 0; i < v.len; i++) {
                norm_squared += v.data[i] * v.data[i];
            }
            for (size_t i = 0; i < P.rows; i++) {
                for (size_t j = 0; j < P.cols; j++) {
                    matrix_set(P, i, j, ((i == j) ? 1.0 : 0.0) - (2.0 * v.data[i] * v.data[j]) / norm_squared);
                }
            }
            // A = P * A * P
            matrix_mul_over(temp, A, P);
            matrix_mul_over(A, P, temp);
            // U = U * P
            matrix_mul_over(temp, U, P);
            matrix_copy_over(U, temp);
        }
    }
    matrix_dealloc(&temp);
    matrix_dealloc(&P);
    vector_dealloc(&v);
}

// Schur decomposition (version with less memory usage)
// In this version, the T matrix is stored in the matrix A, changing its contents!
// The matrix U must be preallocated!
// WARNING! This function changes the contents of A and U!
// Decomposition A = U * T * U^T,
// with T in upper triangular form and U orthogonal
void matrix_schur_dec_over(const Matrix A, const Matrix U) {
    if (!matrix_is_squared(A) || !matrix_is_squared(U) || (A.rows != U.rows)) {
        return;  // Invalid operation
    }
    Matrix Q = matrix_alloc(A.rows, A.cols);
    Matrix R = matrix_alloc(A.rows, A.cols);
    if (matrix_is_valid(Q) && matrix_is_valid(R)) {
        // QR algorithm
        matrix_upper_hessenberg_over(A, U);
        for (size_t i = 0; i < MAX_ITERATIONS; i++) {
            matrix_copy_over(Q, A);
            matrix_qr_dec_over(Q, R);
            // T = R * Q
            matrix_mul_over(A, R, Q);
            // U = U * Q
            matrix_mul_over(R, U, Q);
            matrix_copy_over(U, R);
            if (matrix_is_upper_triangular(A)) {
                break;
            }
        }
    }
    matrix_dealloc(&Q);
    matrix_dealloc(&R);
}

// Matrix diagonalization (version with less memory usage)
// In this version, the D matrix is stored in the matrix A, changing its contents!
// The matrix P must be preallocated!
// WARNING! This function changes the contents of A and P!
// Decomposition A = P * D * P^-1, with D in diagonal form
// Only work for real non repeated eigenvalues
// Based on: https://math.stackexchange.com/questions/3947108/how-to-get-eigenvectors-using-qr-algorithm
void matrix_diagonalization_over(const Matrix A, const Matrix P) {
    if (!matrix_is_squared(A) || !matrix_is_squared(P) || (A.rows != P.rows)) {
        return;  // Invalid operation
    }
    Matrix V = matrix_identity(A.rows);
    Matrix U = matrix_alloc(A.rows, A.cols);
    if (matrix_is_valid(V) && matrix_is_valid(U)) {
        matrix_schur_dec_over(A, U);
        if (matrix_is_upper_triangular(A)) {
            for (size_t i = 1; i < A.rows; i++) {
                for (size_t j = 0; j < i; j++) {
                    matrix_set(V, j, i, matrix_get(A, j, i));
                }
                // Back substitution
                for (size_t j = (i - 1); j < i; j--) {
                    for (size_t k = (j + 1); k < i; k++) {
                        matrix_inc(V, j, i, matrix_get(A, j, k) * matrix_get(V, k, i));
                    }
                    matrix_set(V, j, i, matrix_get(V, j, i) / (matrix_get(A, i, i) - matrix_get(A, j, j)));
                }
            }
            matrix_mul_over(P, U, V);
            // Set the non diagonal elements of A to zero
            for (size_t i = 0; i < A.rows; i++) {
                for (size_t j = (i + 1); j < A.rows; j++) {
                    matrix_set(A, i, j, 0.0);
                    matrix_set(A, j, i, 0.0);
                }
            }
        }
    }
    matrix_dealloc(&U);
    matrix_dealloc(&V);
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