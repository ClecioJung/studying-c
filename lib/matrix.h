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

#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdbool.h>
#include <stdlib.h>

#include "vector.h"

typedef struct {
    double *data;
    size_t rows;
    size_t cols;
} Matrix;

Matrix matrix_alloc(const size_t rows, const size_t cols);
void matrix_dealloc(Matrix *const matrix);
bool matrix_is_valid(const Matrix A);
void matrix_set(const Matrix A, const size_t row, const size_t col, const double value);
double matrix_get(const Matrix A, const size_t row, const size_t col);
void matrix_inc(const Matrix A, const size_t row, const size_t col, const double value);
void matrix_dec(const Matrix A, const size_t row, const size_t col, const double value);
bool matrix_is_squared(const Matrix A);
Matrix matrix_random(const size_t rows, const size_t cols, const double min, const double max);
Matrix matrix_init(const size_t rows, const size_t cols, const double value);
Matrix matrix_identity(const size_t rows);
Matrix matrix_copy(const Matrix matrix);
void matrix_assign(Matrix *const matrix, const Matrix equals);
void matrix_replace(Matrix *const matrix, const Matrix equals);
void matrix_print(const Matrix matrix);
Vector vector_from_matrix_column(const Matrix A, const size_t col);
Matrix matrix_scale(const double scalar, const Matrix A);
Matrix matrix_sum(const Matrix A, const Matrix B);
Matrix matrix_sub(const Matrix A, const Matrix B);
Matrix matrix_mul(const Matrix A, const Matrix B);
Matrix matrix_mul_three(const Matrix A, const Matrix B, const Matrix C);
Vector matrix_mul_vector(const Matrix A, const Vector b);
Matrix matrix_transpose(const Matrix A);
Matrix matrix_symmetric(const Matrix A);
Matrix matrix_skew_symmetric(const Matrix A);
double trace(const Matrix A);
double determinant(const Matrix A);
void lu_decomposition(const Matrix A, Matrix *const L, Matrix *const U);
void lu_crout_decomposition(const Matrix A, Matrix *const L, Matrix *const U);
void qr_decomposition(const Matrix A, Matrix *const Q, Matrix *const R);
Matrix householder_matrix(const Vector vec);
void upper_hessenberg_matrix(const Matrix A, Matrix *const U, Matrix *const H);
void schur_decomposition(const Matrix A, Matrix *const U, Matrix *const T);
Vector eigenvalues(const Matrix A);
double power_method(const Matrix A, Vector *const vec);
Matrix matrix_inverse(const Matrix A);
Matrix pseudo_inverse(const Matrix A);
double matrix_max_diff(const Matrix A, const Matrix previous_A);
bool matrix_are_equal(const Matrix A, const Matrix B);
bool matrix_is_orthogonal(const Matrix A);
bool matrix_is_upper_triangular(const Matrix A);
bool matrix_is_lower_triangular(const Matrix A);
bool matrix_is_symmetric(const Matrix A);
bool matrix_is_skew_symmetric(const Matrix A);
bool matrix_is_null_space(const Matrix A, const Vector vec);

// 'over' methods override the contents of their arguments,
// avoiding the need to allocate more memory for the results
void vector_from_matrix_column_over(const Vector *const vector, const Matrix A, const size_t col);
void matrix_identity_over(const Matrix *const matrix);

#endif  // __MATRIX_H

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