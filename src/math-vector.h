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

#ifndef __MATH_VECTOR_H
#define __MATH_VECTOR_H

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef FUNC_DEF
#define FUNC_DEF static inline
#endif  // FUNC_DEF

#define MAX_ITERATIONS 10000
#define PRECISION 1e-9

typedef double data_type;

typedef struct {
    data_type *data;
    size_t len;
} Vector;

typedef struct {
    data_type *data;
    size_t rows;
    size_t cols;
} Matrix;

// Scalar functions
FUNC_DEF void swap(data_type *const a, data_type *const b);
FUNC_DEF uint64_t factorial(uint64_t value);
FUNC_DEF bool are_close(const data_type a, const data_type b, const data_type delta);
FUNC_DEF data_type maximum(const data_type a, const data_type b);
FUNC_DEF data_type minimum(const data_type a, const data_type b);
FUNC_DEF data_type sign(const data_type value);
FUNC_DEF data_type random_number(const data_type min, const data_type max);
FUNC_DEF data_type square_root(const data_type value);
FUNC_DEF data_type power(const data_type base, uint64_t expoent);
FUNC_DEF data_type exponential(const data_type value);

// Vector functions
FUNC_DEF Vector vector_alloc(const size_t len);
FUNC_DEF void vector_dealloc(Vector *const vector);
FUNC_DEF Vector vector_realloc(Vector *const vector, const size_t len);
FUNC_DEF bool vector_is_valid(const Vector vector);
FUNC_DEF void vector_set(const Vector vector, const size_t index, const data_type value);
FUNC_DEF data_type vector_get(const Vector vector, const size_t index);
FUNC_DEF Vector vector_random(const size_t len, const data_type min, const data_type max);
FUNC_DEF Vector vector_init(const size_t len, const data_type value);
FUNC_DEF void vector_assign(Vector *const vector, const Vector equals);
FUNC_DEF void vector_replace(Vector *const vector, const Vector equals);
FUNC_DEF Vector vector_copy(const Vector vector);
FUNC_DEF void vector_print(const Vector vector);
FUNC_DEF Vector vector_scale(const data_type scalar, const Vector vector);
FUNC_DEF data_type dot_product(const Vector a, const Vector b);
FUNC_DEF Vector cross_product(const Vector a, const Vector b);
FUNC_DEF data_type euclidean_norm(const Vector x);
FUNC_DEF data_type vector_max(const Vector x);
FUNC_DEF Vector vector_sum(const Vector a, const Vector b);
FUNC_DEF Vector vector_sub(const Vector a, const Vector b);
FUNC_DEF data_type vector_max_diff(Vector x, Vector previous_x);
FUNC_DEF bool vector_are_equal(const Vector a, const Vector b);

// Sorting functions
FUNC_DEF bool vector_is_sorted(const Vector vec);
FUNC_DEF void bubble_sort(const Vector vec);
FUNC_DEF void select_sort(const Vector vec);
FUNC_DEF void insert_sort(const Vector vec);
FUNC_DEF void shell_sort(const Vector vec);
FUNC_DEF void heap_sort(const Vector vec);
FUNC_DEF void merge_sort(const Vector vec);
FUNC_DEF void quicksort(const Vector vec);

// Search functions
FUNC_DEF size_t sequential_search(const Vector vec, const data_type value);
FUNC_DEF size_t binary_search(const Vector vec, const data_type value);

// Matrix functions
FUNC_DEF Matrix matrix_alloc(const size_t rows, const size_t cols);
FUNC_DEF void matrix_dealloc(Matrix *const matrix);
FUNC_DEF bool matrix_is_valid(const Matrix A);
FUNC_DEF void matrix_set(const Matrix A, const size_t row, const size_t col, const data_type value);
FUNC_DEF data_type matrix_get(const Matrix A, const size_t row, const size_t col);
FUNC_DEF bool matrix_is_squared(const Matrix A);
FUNC_DEF Matrix matrix_random(const size_t rows, const size_t cols, const data_type min, const data_type max);
FUNC_DEF Matrix matrix_init(const size_t rows, const size_t cols, const data_type value);
FUNC_DEF Matrix matrix_identity(const size_t rows);
FUNC_DEF Matrix matrix_copy(const Matrix matrix);
FUNC_DEF void matrix_assign(Matrix *const matrix, const Matrix equals);
FUNC_DEF void matrix_replace(Matrix *const matrix, const Matrix equals);
FUNC_DEF void matrix_print(const Matrix matrix);
FUNC_DEF Vector vector_from_matrix_column(const Matrix A, const size_t col);
FUNC_DEF Matrix matrix_scale(const data_type scalar, const Matrix A);
FUNC_DEF Matrix matrix_mul(const Matrix A, const Matrix B);
FUNC_DEF Matrix matrix_mul_three(const Matrix A, const Matrix B, const Matrix C);
FUNC_DEF Vector matrix_mul_vector(const Matrix A, const Vector b);
FUNC_DEF Matrix matrix_transpose(const Matrix A);
FUNC_DEF Matrix matrix_symmetric(const Matrix A);
FUNC_DEF Matrix matrix_skew_symmetric(const Matrix A);
FUNC_DEF data_type trace(const Matrix A);
FUNC_DEF data_type determinant(const Matrix A);
FUNC_DEF void lu_decomposition(const Matrix A, Matrix *const L, Matrix *const U);
FUNC_DEF void lu_crout_decomposition(const Matrix A, Matrix *const L, Matrix *const U);
FUNC_DEF void qr_decomposition(const Matrix A, Matrix *const Q, Matrix *const R);
FUNC_DEF Matrix householder_matrix(const Vector vec);
FUNC_DEF void upper_hessenberg_matrix(const Matrix A, Matrix *const U, Matrix *const H);
FUNC_DEF void schur_decomposition(const Matrix A, Matrix *const U, Matrix *const T);
FUNC_DEF Vector eigenvalues(const Matrix A);
FUNC_DEF data_type power_method(const Matrix A, Vector *const vec);
FUNC_DEF Matrix matrix_inverse(const Matrix A);
FUNC_DEF Matrix pseudo_inverse(const Matrix A);
FUNC_DEF Matrix matrix_sum(const Matrix A, const Matrix B);
FUNC_DEF Matrix matrix_sub(const Matrix A, const Matrix B);
FUNC_DEF data_type matrix_max_diff(const Matrix A, const Matrix previous_A);
FUNC_DEF bool matrix_are_equal(const Matrix A, const Matrix B);
FUNC_DEF bool matrix_is_orthogonal(const Matrix A);
FUNC_DEF bool matrix_is_upper_triangular(const Matrix A);

// Methods for solving linear systems
FUNC_DEF Vector back_substitution(const Matrix A, const Vector b);
FUNC_DEF Vector forward_substitution(const Matrix A, const Vector b);
FUNC_DEF Vector gaussian_elimination(const Matrix A, const Vector b);
FUNC_DEF Vector gauss_jordan(const Matrix A, const Vector b);
FUNC_DEF Vector lu_solving(const Matrix A, const Vector b);
FUNC_DEF Vector jacobi_method(const Matrix A, const Vector b);
FUNC_DEF Vector gauss_seidel(const Matrix A, const Vector b);
FUNC_DEF bool columns_condition(const Matrix A);
FUNC_DEF bool rows_condition(const Matrix A);
FUNC_DEF bool sassenfeld_condition(const Matrix A);

// Curve fitting functions
FUNC_DEF data_type lagrange_interpolation(const Vector x, const Vector y, const data_type value);
FUNC_DEF data_type linear_regression(const Vector x, const Vector y, data_type *const a, data_type *const b);
FUNC_DEF Vector polynomial_regression(const Vector x, const Vector y, const size_t order);
FUNC_DEF data_type compute_polynomial(const Vector coefficients, const data_type x);

#endif  // __MATH_VECTOR_H

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------

#ifdef MATH_VECTOR_IMPLEMENTATION

FUNC_DEF void swap(data_type *const a, data_type *const b) {
    if ((a != NULL) && (b != NULL)) {
        const data_type temp = *a;
        *a = *b;
        *b = temp;
    }
}

FUNC_DEF uint64_t factorial(uint64_t value) {
    uint64_t result = 1;
    while (value > 1) {
        result = result * value;
        value--;
    }
    return result;
}

FUNC_DEF bool are_close(const data_type a, const data_type b, const data_type delta) {
    return (fabs(a - b) < delta);
}

FUNC_DEF data_type maximum(const data_type a, const data_type b) {
    return ((a > b) ? a : b);
}

FUNC_DEF data_type minimum(const data_type a, const data_type b) {
    return ((a < b) ? a : b);
}

FUNC_DEF data_type sign(const data_type value) {
    return ((value > 0) ? 1.0 : -1.0);
}

FUNC_DEF data_type random_number(const data_type min, const data_type max) {
    const data_type rand_unitary = ((data_type)rand()) / ((data_type)RAND_MAX);
    return rand_unitary * (max - min) + min;
}

// My own square root function, so I don't need to link with -lm,
// avoiding any dependencies
FUNC_DEF data_type square_root(const data_type value) {
    if (fabs(value) < PRECISION) {
        return 0.0;  // square_root(0) = 0
    } else if (value < 0) {
        return NAN;
    }
    data_type x = value;
    // Use the Newton-Raphson method to find the root
    // of the function f(x) = x^2 - value
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const data_type delta = (value / x - x) / 2.0;
        x += delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return x;
}

// My own power function, so I don't need to link with -lm,
// avoiding any dependencies
FUNC_DEF data_type power(const data_type base, uint64_t expoent) {
    data_type result = 1.0;
    while (expoent != 0) {
        result *= base;
        expoent--;
    }
    return result;
}

// My own exponential function, so I don't need to link with -lm,
// avoiding any dependencies
FUNC_DEF data_type exponential(const data_type value) {
    data_type term = 1.0;
    data_type result = 1.0;
    for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
        term *= value / ((data_type)i);
        result += term;
        if (fabs(term) < PRECISION) {
            break;
        }
    }
    return result;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector vector_alloc(const size_t len) {
    if (len == 0) {
        return (Vector){0};
    }
    Vector vector = (Vector){
        .data = (data_type *)malloc(len * sizeof(data_type)),
        .len = len,
    };
    if (vector.data == NULL) {
        vector.len = 0;
    }
    return vector;
}

FUNC_DEF void vector_dealloc(Vector *const vector) {
    if (vector == NULL) {
        return;
    }
    free(vector->data);
    vector->data = NULL;
    vector->len = 0;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector vector_realloc(Vector *const vector, const size_t len) {
    if (len != 0) {
        data_type *new_data = (data_type *)realloc(vector->data, len * sizeof(data_type));
        if (new_data != NULL) {
            vector->len = len;
            vector->data = new_data;
        }
    } else {
        vector_dealloc(vector);
    }
    return *vector;
}

FUNC_DEF bool vector_is_valid(const Vector vector) {
    return (vector.data != NULL);
}

FUNC_DEF void vector_set(const Vector vector, const size_t index, const data_type value) {
    if (vector_is_valid(vector) || (index >= vector.len)) {
        return;
    }
    vector.data[index] = value;
}

FUNC_DEF data_type vector_get(const Vector vector, const size_t index) {
    if (vector_is_valid(vector) || (index >= vector.len)) {
        return NAN;
    }
    return vector.data[index];
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector vector_random(const size_t len, const data_type min, const data_type max) {
    Vector vector = vector_alloc(len);
    if (vector_is_valid(vector)) {
        srand(time(NULL));
        for (size_t i = 0; i < len; i++) {
            vector.data[i] = random_number(min, max);
        }
    }
    return vector;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector vector_init(const size_t len, const data_type value) {
    Vector vector = vector_alloc(len);
    if (vector_is_valid(vector)) {
        for (size_t i = 0; i < len; i++) {
            vector.data[i] = value;
        }
    }
    return vector;
}

FUNC_DEF void vector_assign(Vector *const vector, const Vector equals) {
    vector->data = equals.data;
    vector->len = equals.len;
}

FUNC_DEF void vector_replace(Vector *const vector, const Vector equals) {
    vector_dealloc(vector);
    vector_assign(vector, equals);
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector vector_copy(const Vector vector) {
    Vector new_vec = vector_alloc(vector.len);
    if (vector_is_valid(new_vec)) {
        for (size_t i = 0; i < vector.len; i++) {
            new_vec.data[i] = vector.data[i];
        }
    }
    return new_vec;
}

FUNC_DEF void vector_print(const Vector vector) {
    for (size_t i = 0; i < vector.len; i++) {
        const data_type value = fabs(vector.data[i]) > PRECISION ? vector.data[i] : 0.0;
        printf("[%zu]: %g\n", i, value);
    }
    printf("\n");
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector vector_scale(const data_type scalar, const Vector vector) {
    Vector new_vec = vector_alloc(vector.len);
    if (vector_is_valid(new_vec)) {
        for (size_t i = 0; i < vector.len; i++) {
            new_vec.data[i] = scalar * vector.data[i];
        }
    }
    return new_vec;
}

FUNC_DEF data_type dot_product(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return 0.0;  // Invalid operation
    }
    data_type result = 0.0;
    for (size_t i = 0; i < a.len; i++) {
        result += a.data[i] * b.data[i];
    }
    return result;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector cross_product(const Vector a, const Vector b) {
    if ((a.len != 3) || (b.len != 3)) {
        return (Vector){0};  // Invalid operation
    }
    Vector result = vector_alloc(a.len);
    result.data[0] = a.data[1] * b.data[2] - a.data[2] * b.data[1];
    result.data[1] = a.data[2] * b.data[0] - a.data[0] * b.data[2];
    result.data[2] = a.data[0] * b.data[1] - a.data[1] * b.data[0];
    return result;
}

FUNC_DEF data_type euclidean_norm(const Vector x) {
    data_type value = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        value += x.data[i] * x.data[i];
    }
    return square_root(value);
}

FUNC_DEF data_type vector_max(const Vector x) {
    data_type value = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        value = maximum(value, x.data[i]);
    }
    return value;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector vector_sum(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector new_vec = vector_alloc(a.len);
    if (vector_is_valid(new_vec)) {
        for (size_t i = 0; i < a.len; i++) {
            new_vec.data[i] = a.data[i] + b.data[i];
        }
    }
    return new_vec;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector vector_sub(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector new_vec = vector_alloc(a.len);
    if (vector_is_valid(new_vec)) {
        for (size_t i = 0; i < a.len; i++) {
            new_vec.data[i] = a.data[i] - b.data[i];
        }
    }
    return new_vec;
}

FUNC_DEF data_type vector_max_diff(const Vector x, const Vector previous_x) {
    if (x.len != previous_x.len) {
        return INFINITY;  // Invalid operation
    }
    data_type error = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        error = maximum(fabs(x.data[i] - previous_x.data[i]), error);
    }
    return error;
}

FUNC_DEF bool vector_are_equal(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return false;
    }
    for (size_t i = 0; i < a.len; i++) {
        if (!are_close(a.data[i], b.data[i], PRECISION)) {
            return false;
        }
    }
    return true;
}

FUNC_DEF bool vector_is_sorted(const Vector vec) {
    for (size_t i = 0; (i + 1) < vec.len; i++) {
        if (vec.data[i] > vec.data[i + 1]) {
            return false;
        }
    }
    return true;
}

FUNC_DEF void bubble_sort(const Vector vec) {
    for (size_t i = 1; i < vec.len; i++) {
        for (size_t j = (vec.len - 1); j >= i; j--) {
            if (vec.data[j - 1] > vec.data[j]) {
                swap(&vec.data[j - 1], &vec.data[j]);
            }
        }
    }
}

FUNC_DEF void select_sort(const Vector vec) {
    for (size_t i = 0; (i + 1) < vec.len; i++) {
        size_t k = i;
        data_type temp = vec.data[i];
        for (size_t j = (i + 1); j < vec.len; j++) {
            if (vec.data[j] < temp) {
                k = j;
                temp = vec.data[j];
            }
        }
        if (k != i) {
            vec.data[k] = vec.data[i];
            vec.data[i] = temp;
        }
    }
}

FUNC_DEF void insert_sort(const Vector vec) {
    for (size_t i = 1; i < vec.len; i++) {
        const data_type temp = vec.data[i];
        size_t j = (i - 1);
        while (j < vec.len && temp < vec.data[j]) {
            vec.data[j + 1] = vec.data[j];
            j--;
        }
        vec.data[j + 1] = temp;
    }
}

FUNC_DEF void shell_sort(const Vector vec) {
    // Ciura gap sequence
    static const size_t gaps[] = {701, 301, 132, 57, 23, 10, 4, 1};
    static const size_t gap_len = sizeof(gaps) / sizeof(gaps[0]);

    for (size_t k = 0; k < gap_len; k++) {
        const size_t gap = gaps[k];
        for (size_t i = gap; i < vec.len; i++) {
            const data_type temp = vec.data[i];
            size_t j = (i - gap);
            while (j < vec.len && temp < vec.data[j]) {
                vec.data[j + gap] = vec.data[j];
                j = (j - gap);
            }
            vec.data[j + gap] = temp;
        }
    }
}

// Heapify a subtree rooted with node in array at index
FUNC_DEF void heapify(data_type *const array, const size_t len, const size_t index) {
    size_t largest = index;
    const size_t left = 2 * index + 1;
    const size_t right = 2 * index + 2;
    // If left child is larger than root
    if (left < len && array[left] > array[largest]) {
        largest = left;
    }
    // If right child is larger than largest so far
    if (right < len && array[right] > array[largest]) {
        largest = right;
    }
    // If largest is not root
    if (largest != index) {
        swap(&array[largest], &array[index]);
        // Recursively heapify the affected sub-tree
        heapify(array, len, largest);
    }
}
FUNC_DEF void heap_sort(const Vector vec) {
    // Build heap (rearrange array)
    for (size_t i = vec.len / 2 - 1; i < vec.len; i--) {
        heapify(vec.data, vec.len, i);
    }
    // One by one extract an element from heap
    for (size_t i = vec.len - 1; i > 0; i--) {
        swap(&vec.data[0], &vec.data[i]);
        // Call max heapify on the reduced heap
        heapify(vec.data, i, 0);
    }
}

// Merge two subarrays L and M into A
FUNC_DEF void merge(data_type *const array, const size_t p, const size_t q, const size_t r) {
    // Create L <- A[p..q] and M <- A[q+1..r]
    const size_t n1 = q - p + 1;
    const size_t n2 = r - q;
    data_type L[n1], M[n2];
    for (size_t i = 0; i < n1; i++) {
        L[i] = array[p + i];
    }
    for (size_t j = 0; j < n2; j++) {
        M[j] = array[q + 1 + j];
    }
    // Maintain current index of sub-arrays and main array
    size_t i = 0;
    size_t j = 0;
    size_t k = p;
    // Until we reach either end of either L or M, pick larger among
    // elements L and M and place them in the correct position at A[p..r]
    while (i < n1 && j < n2) {
        if (L[i] <= M[j]) {
            array[k] = L[i];
            i++;
        } else {
            array[k] = M[j];
            j++;
        }
        k++;
    }
    // When we run out of elements in either L or M,
    // pick up the remaining elements and put in A[p..r]
    while (i < n1) {
        array[k] = L[i];
        i++;
        k++;
    }
    while (j < n2) {
        array[k] = M[j];
        j++;
        k++;
    }
}
// Divide the array into two subarrays, sort them and merge them
FUNC_DEF void split_and_sort(data_type *const array, const size_t left, const size_t right) {
    if (left < right) {
        data_type middle = left + (right - left) / 2;
        split_and_sort(array, left, middle);
        split_and_sort(array, middle + 1, right);
        merge(array, left, middle, right);
    }
}
FUNC_DEF void merge_sort(const Vector vec) {
    if (vec.len > 0) {
        split_and_sort(vec.data, 0, vec.len - 1);
    }
}

FUNC_DEF void qs(data_type *const array, const size_t left, const size_t right) {
    const data_type middle = array[(left + right) / 2];
    size_t i = left;
    size_t j = right;
    do {
        while (array[i] < middle && i < right) {
            i++;
        }
        while (middle < array[j] && j > left && j != 0) {
            j--;
        }
        if (i <= j) {
            swap(&array[i], &array[j]);
            i++;
            if (j != 0) {
                j--;
            }
        }
    } while (i <= j);
    if (left < j) {
        qs(array, left, j);
    }
    if (i < right) {
        qs(array, i, right);
    }
}
FUNC_DEF void quicksort(const Vector vec) {
    if (vec.len > 0) {
        qs(vec.data, 0, vec.len - 1);
    }
}

FUNC_DEF size_t sequential_search(const Vector vec, const data_type value) {
    size_t closest_index = 0;
    for (size_t i = 0; i < vec.len; i++) {
        if (vec.data[i] == value) {
            return i;
        } else if (fabs(vec.data[i] - value) < fabs(vec.data[closest_index] - value)) {
            closest_index = i;
        }
    }
    return closest_index;
}

// It works only with ordered arrays
FUNC_DEF size_t binary_search(const Vector vec, const data_type value) {
    size_t middle = vec.len;
    size_t low = 0;
    size_t high = vec.len - 1;
    while (low <= high) {
        middle = (low + high) / 2;
        if (value < vec.data[middle]) {
            high = middle - 1;
        } else if (value > vec.data[middle]) {
            low = middle + 1;
        } else {
            break;
        }
    }
    return middle;
}

// Remember to free the matrix after calling this function!
FUNC_DEF Matrix matrix_alloc(const size_t rows, const size_t cols) {
    if ((rows == 0) || (cols == 0)) {
        return (Matrix){0};
    }
    Matrix matrix = (Matrix){
        .data = (data_type *)malloc(rows * cols * sizeof(data_type)),
        .rows = rows,
        .cols = cols,
    };
    if (matrix.data == NULL) {
        matrix.rows = 0;
        matrix.cols = 0;
    }
    return matrix;
}

FUNC_DEF void matrix_dealloc(Matrix *const matrix) {
    if (matrix == NULL) {
        return;
    }
    free(matrix->data);
    matrix->data = NULL;
    matrix->rows = 0;
    matrix->cols = 0;
}

FUNC_DEF bool matrix_is_valid(const Matrix A) {
    return (A.data != NULL);
}

FUNC_DEF void matrix_set(const Matrix A, const size_t row, const size_t col, const data_type value) {
    if (!matrix_is_valid(A) || (row >= A.rows) || (col >= A.cols)) {
        return;
    }
    A.data[row * A.cols + col] = value;
}

FUNC_DEF data_type matrix_get(const Matrix A, const size_t row, const size_t col) {
    if (!matrix_is_valid(A) || (row >= A.rows) || (col >= A.cols)) {
        return NAN;
    }
    return A.data[row * A.cols + col];
}

FUNC_DEF bool matrix_is_squared(const Matrix A) {
    return (matrix_is_valid(A) && (A.rows == A.cols));
}

// Remember to free the matrix after calling this function!
FUNC_DEF Matrix matrix_random(const size_t rows, const size_t cols, const data_type min, const data_type max) {
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
FUNC_DEF Matrix matrix_init(const size_t rows, const size_t cols, const data_type value) {
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
FUNC_DEF Matrix matrix_identity(const size_t rows) {
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
FUNC_DEF Matrix matrix_copy(const Matrix matrix) {
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

FUNC_DEF void matrix_assign(Matrix *const matrix, const Matrix equals) {
    matrix->data = equals.data;
    matrix->rows = equals.rows;
    matrix->cols = equals.cols;
}

FUNC_DEF void matrix_replace(Matrix *const matrix, const Matrix equals) {
    matrix_dealloc(matrix);
    matrix_assign(matrix, equals);
}

FUNC_DEF void matrix_print(const Matrix matrix) {
    printf("%*s", 2, "");
    for (size_t j = 0; j < matrix.cols; j++) {
        printf("%*s[%zu] ", 5, "", j);
    }
    printf("\n");
    for (size_t i = 0; i < matrix.rows; i++) {
        printf("[%zu]: ", i);
        for (size_t j = 0; j < matrix.cols; j++) {
            const data_type value = fabs(matrix_get(matrix, i, j)) > PRECISION ? matrix_get(matrix, i, j) : 0.0;
            printf("%-10g ", value);
        }
        printf("\n");
    }
    printf("\n");
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector vector_from_matrix_column(const Matrix A, const size_t col) {
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
FUNC_DEF Matrix matrix_scale(const data_type scalar, const Matrix A) {
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
FUNC_DEF Matrix matrix_mul(const Matrix A, const Matrix B) {
    if (A.cols != B.rows) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = matrix_init(A.rows, B.cols, 0.0);
    if (matrix_is_valid(result)) {
        for (size_t i = 0; i < result.rows; i++) {
            for (size_t j = 0; j < result.cols; j++) {
                for (size_t k = 0; k < A.cols; k++) {
                    matrix_set(result, i, j, matrix_get(result, i, j) + matrix_get(A, i, k) * matrix_get(B, k, j));
                }
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
FUNC_DEF Matrix matrix_mul_three(const Matrix A, const Matrix B, const Matrix C) {
    if ((A.cols != B.rows) || (B.cols != C.rows)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = matrix_mul(A, B);
    matrix_replace(&result, matrix_mul(result, C));
    return result;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector matrix_mul_vector(const Matrix A, const Vector b) {
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
FUNC_DEF Matrix matrix_transpose(const Matrix A) {
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
FUNC_DEF Matrix matrix_symmetric(const Matrix A) {
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
FUNC_DEF Matrix matrix_skew_symmetric(const Matrix A) {
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

FUNC_DEF data_type trace(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return 0.0;  // Invalid operation
    }
    data_type trace = 0.0;
    for (size_t i = 0; i < A.rows; i++) {
        trace += matrix_get(A, i, i);
    }
    return trace;
}

FUNC_DEF data_type determinant(const Matrix A) {
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
            const data_type m = matrix_get(B, i, k) / matrix_get(B, k, k);
            for (size_t j = (k + 1); j < B.rows; j++) {
                matrix_set(B, i, j, matrix_get(B, i, j) - m * matrix_get(B, k, j));
            }
        }
    }
    data_type det = 1.0;
    for (size_t k = 0; k < B.rows; k++) {
        det *= matrix_get(B, k, k);
    }
    matrix_dealloc(&B);
    return det;
}

// Remember to free L and U matrices after calling this function!
FUNC_DEF void lu_decomposition(const Matrix A, Matrix *const L, Matrix *const U) {
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
                matrix_set(*U, i, j, matrix_get(*U, i, j) - matrix_get(*L, i, k) * matrix_get(*U, k, j));
            }
        }
    }
}

// Remember to free L and U matrices after calling this function!
FUNC_DEF void lu_crout_decomposition(const Matrix A, Matrix *const L, Matrix *const U) {
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
                matrix_set(*L, i, k, matrix_get(*L, i, k) - matrix_get(*L, i, l) * matrix_get(*U, l, k));
            }
        }
        if ((k + 1) != A.rows) {
            for (size_t j = (k + 1); j < A.rows; j++) {
                matrix_set(*U, k, j, matrix_get(A, k, j));
                for (size_t l = 0; l < k; l++) {
                    matrix_set(*U, k, j, matrix_get(*U, k, j) - matrix_get(*L, k, l) * matrix_get(*U, l, j));
                }
                matrix_set(*U, k, j, matrix_get(*U, k, j) / matrix_get(*L, k, k));
            }
        }
    }
}

// Remember to free Q and R matrices after calling this function!
FUNC_DEF void qr_decomposition(const Matrix A, Matrix *const Q, Matrix *const R) {
    if (!matrix_is_squared(A) || (Q == NULL) || (R == NULL)) {
        return;  // Invalid operation
    }
    *Q = matrix_copy(A);
    *R = matrix_init(A.rows, A.cols, 0.0);
    if (!matrix_is_valid(*Q) || !matrix_is_valid(*R)) {
    qr_decomposition_safe_exit:
        matrix_dealloc(Q);
        matrix_dealloc(R);
        return;
    }
    for (size_t i = 0; i < A.cols; i++) {
        Vector ai = vector_from_matrix_column(A, i);
        if (!vector_is_valid(ai)) {
            goto qr_decomposition_safe_exit;
        }
        // Gram-Schmidt process
        for (size_t j = 0; j < i; j++) {
            Vector qj = vector_from_matrix_column(*Q, j);
            if (!vector_is_valid(qj)) {
                vector_dealloc(&ai);
                goto qr_decomposition_safe_exit;
            }
            matrix_set(*R, j, i, dot_product(qj, ai));
            for (size_t k = 0; k < A.rows; k++) {
                matrix_set(*Q, k, i, matrix_get(*Q, k, i) - matrix_get(*R, j, i) * qj.data[k]);
            }
            vector_dealloc(&qj);
        }
        vector_dealloc(&ai);
        Vector qi = vector_from_matrix_column(*Q, i);
        if (!vector_is_valid(qi)) {
            goto qr_decomposition_safe_exit;
        }
        matrix_set(*R, i, i, euclidean_norm(qi));
        for (size_t k = 0; k < A.rows; k++) {
            matrix_set(*Q, k, i, qi.data[k] / matrix_get(*R, i, i));
        }
        vector_dealloc(&qi);
    }
}

// Remember to free the returned matrix after calling this function!
FUNC_DEF Matrix householder_matrix(const Vector vec) {
    Matrix matrix = matrix_identity(vec.len);
    data_type norm_squared = 0.0;
    for (size_t i = 0; i < vec.len; i++) {
        norm_squared += vec.data[i] * vec.data[i];
    }
    for (size_t i = 0; i < matrix.rows; i++) {
        for (size_t j = 0; j < matrix.cols; j++) {
            matrix_set(matrix, i, j, matrix_get(matrix, i, j) - (2.0 * vec.data[i] * vec.data[j]) / norm_squared);
        }
    }
    return matrix;
}

// Remember to free the matrices U and H after calling this function!
// Decomposition A = U * H * U^T,
// with H in upper Hessenberg form and U orthogonal
FUNC_DEF void upper_hessenberg_matrix(const Matrix A, Matrix *const U, Matrix *const H) {
    if (!matrix_is_squared(A) || (U == NULL) || (H == NULL)) {
        return;  // Invalid operation
    }
    *U = matrix_identity(A.rows);
    *H = matrix_copy(A);
    if (!matrix_is_valid(*H) || !matrix_is_valid(*U)) {
    hessenberg_matrix_safe_exit:
        matrix_dealloc(H);
        matrix_dealloc(U);
        return;
    }
    for (size_t i = 1; (i + 1) < A.rows; i++) {
        Vector v = vector_alloc(H->rows - i);
        for (size_t j = 0; j < v.len; j++) {
            v.data[j] = matrix_get(*H, i + j, i - 1);
        }
        v.data[0] += sign(v.data[0]) * euclidean_norm(v);
        Matrix householder = householder_matrix(v);
        vector_dealloc(&v);
        if (!matrix_is_valid(householder)) {
            goto hessenberg_matrix_safe_exit;
        }
        Matrix P = matrix_identity(H->rows);
        if (!matrix_is_valid(P)) {
            matrix_dealloc(&householder);
            goto hessenberg_matrix_safe_exit;
        }
        const size_t offset = H->rows - householder.rows;
        for (size_t i = offset; i < P.rows; i++) {
            for (size_t j = offset; j < P.cols; j++) {
                matrix_set(P, i, j, matrix_get(householder, i - offset, j - offset));
            }
        }
        matrix_dealloc(&householder);
        matrix_replace(H, matrix_mul_three(P, *H, P));
        matrix_replace(U, matrix_mul(*U, P));
        matrix_dealloc(&P);
        if (!matrix_is_valid(*H) || !matrix_is_valid(*U)) {
            goto hessenberg_matrix_safe_exit;
        }
    }
}

// Remember to free the matrices U and T after calling this function!
FUNC_DEF void schur_decomposition(const Matrix A, Matrix *const U, Matrix *const T) {
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
FUNC_DEF data_type power_method(const Matrix A, Vector *const vec) {
    if (!matrix_is_squared(A) || (vec == NULL)) {
        return NAN;  // Invalid operation
    }
    *vec = vector_random(A.rows, 0.0, 1.0);
    if (!vector_is_valid(*vec)) {
    power_method_safe_exit:
        vector_dealloc(vec);
        return NAN;
    }
    data_type eig = 0.0;
    for (size_t i = 0; i < MAX_ITERATIONS; i++) {
        Vector previous_vec;
        vector_assign(&previous_vec, *vec);
        *vec = matrix_mul_vector(A, *vec);
        eig = euclidean_norm(*vec);
        vector_replace(vec, vector_scale((1.0 / eig), *vec));
        if (!vector_is_valid(*vec) || !vector_is_valid(previous_vec)) {
            vector_dealloc(&previous_vec);
            goto power_method_safe_exit;
        }
        const data_type error = vector_max_diff(*vec, previous_vec);
        vector_dealloc(&previous_vec);
        if (error < PRECISION) {
            break;
        }
    }
    return eig;
}

// Remember to free the returned matrix after calling this function!
FUNC_DEF Matrix matrix_inverse(const Matrix A) {
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

FUNC_DEF Matrix pseudo_inverse(const Matrix A) {
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

// Remember to free the returned matrix after calling this function!
FUNC_DEF Matrix matrix_sum(const Matrix A, const Matrix B) {
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
FUNC_DEF Matrix matrix_sub(const Matrix A, const Matrix B) {
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

FUNC_DEF data_type matrix_max_diff(const Matrix A, const Matrix previous_A) {
    if ((A.rows != previous_A.rows) || (A.cols != previous_A.cols)) {
        return INFINITY;  // Invalid operation
    }
    data_type error = 0.0;
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            error = maximum(fabs(matrix_get(A, i, j) - matrix_get(previous_A, i, j)), error);
        }
    }
    return error;
}

FUNC_DEF bool matrix_are_equal(const Matrix A, const Matrix B) {
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

FUNC_DEF bool matrix_is_orthogonal(const Matrix A) {
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

FUNC_DEF bool matrix_is_upper_triangular(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    data_type error = 0.0;
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < i; j++) {
            error = maximum(error, fabs(matrix_get(A, i, j)));
        }
    }
    return (error < PRECISION);
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector back_substitution(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = vector_alloc(b.len);
    if (vector_is_valid(x)) {
        for (size_t i = (b.len - 1); i < b.len; i--) {
            data_type sum = b.data[i];
            for (size_t j = (i + 1); j < b.len; j++) {
                sum -= matrix_get(A, i, j) * x.data[j];
            }
            x.data[i] = sum / matrix_get(A, i, i);
        }
    }
    return x;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector forward_substitution(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = vector_alloc(b.len);
    if (vector_is_valid(x)) {
        for (size_t i = 0; i < b.len; i++) {
            data_type sum = b.data[i];
            for (size_t j = 0; j < i; j++) {
                sum -= matrix_get(A, i, j) * x.data[j];
            }
            x.data[i] = sum / matrix_get(A, i, i);
        }
    }
    return x;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector gaussian_elimination(const Matrix A, const Vector b) {
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
        data_type w = fabs(matrix_get(A_copied, k, k));
        size_t r = k;
        for (size_t i = (k + 1); i < A_copied.rows; i++) {
            if (fabs(matrix_get(A_copied, i, k)) > w) {
                w = fabs(matrix_get(A_copied, i, k));
                r = i;
            }
        }
        if (r != k) {
            for (size_t i = k; i < A_copied.rows; i++) {
                const data_type temp = matrix_get(A_copied, k, i);
                matrix_set(A_copied, k, i, matrix_get(A_copied, r, i));
                matrix_set(A_copied, r, i, temp);
            }
            swap(&b_copied.data[k], &b_copied.data[r]);
        }
    }
    // Gaussian elimination
    for (size_t k = 0; (k + 1) < A_copied.rows; k++) {
        for (size_t i = (k + 1); i < A_copied.rows; i++) {
            const data_type m = matrix_get(A_copied, i, k) / matrix_get(A_copied, k, k);
            for (size_t j = (k + 1); j < A_copied.rows; j++) {
                matrix_set(A_copied, i, j, matrix_get(A_copied, i, j) - m * matrix_get(A_copied, k, j));
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
FUNC_DEF Vector gauss_jordan(const Matrix A, const Vector b) {
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
            const data_type m = matrix_get(A_copied, i, k) / matrix_get(A_copied, k, k);
            const data_type p = matrix_get(A_copied, k, k);
            if (i == k) {
                for (size_t j = k; j < A_copied.cols; j++) {
                    matrix_set(A_copied, i, j, matrix_get(A_copied, i, j) / p);
                }
                b_copied.data[i] /= p;
            } else {
                for (size_t j = k; j < A_copied.cols; j++) {
                    matrix_set(A_copied, i, j, matrix_get(A_copied, i, j) - m * matrix_get(A_copied, k, j));
                }
                b_copied.data[i] -= m * b_copied.data[k];
            }
        }
    }
    matrix_dealloc(&A_copied);
    return b_copied;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector lu_solving(const Matrix A, const Vector b) {
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
FUNC_DEF Vector jacobi_method(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = vector_init(b.len, 0.0);
    if (vector_is_valid(x)) {
        for (size_t k = 0; k < MAX_ITERATIONS; k++) {
            Vector previous_x = vector_copy(x);
            if (!vector_is_valid(previous_x)) {
                vector_dealloc(&x);
                break;
            }
            for (size_t i = 0; i < x.len; i++) {
                data_type sum = b.data[i];
                for (size_t j = 0; j < x.len; j++) {
                    if (i != j) {
                        sum -= matrix_get(A, i, j) * previous_x.data[j];
                    }
                }
                x.data[i] = sum / matrix_get(A, i, i);
            }
            const data_type error = vector_max_diff(x, previous_x);
            vector_dealloc(&previous_x);
            if (error < PRECISION) {
                break;
            }
        }
    }
    return x;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector gauss_seidel(const Matrix A, const Vector b) {
    if (!matrix_is_squared(A) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = vector_init(b.len, 0.0);
    if (vector_is_valid(x)) {
        for (size_t k = 0; k < MAX_ITERATIONS; k++) {
            Vector previous_x = vector_copy(x);
            if (!vector_is_valid(previous_x)) {
                vector_dealloc(&x);
                break;
            }
            for (size_t i = 0; i < b.len; i++) {
                data_type sum = b.data[i];
                for (size_t j = 0; j < x.len; j++) {
                    if (i != j) {
                        sum -= matrix_get(A, i, j) * x.data[j];
                    }
                }
                x.data[i] = sum / matrix_get(A, i, i);
            }
            const data_type error = vector_max_diff(x, previous_x);
            vector_dealloc(&previous_x);
            if (error < PRECISION) {
                break;
            }
        }
    }
    return x;
}

FUNC_DEF bool columns_condition(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        data_type sum = 0.0;
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

FUNC_DEF bool rows_condition(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        data_type sum = 0.0;
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

FUNC_DEF bool sassenfeld_condition(const Matrix A) {
    if (!matrix_is_squared(A)) {
        return false;  // Invalid operation
    }
    bool condition = false;
    Vector beta = vector_alloc(A.rows);
    if (vector_is_valid(beta)) {
        for (size_t i = 0; i < A.rows; i++) {
            data_type sum = 0.0;
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

FUNC_DEF data_type lagrange_interpolation(const Vector x, const Vector y, const data_type value) {
    if (x.len != y.len) {
        return 0.0;  // Invalid operation
    }
    data_type sum = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        data_type product = y.data[i];
        for (size_t j = 0; j < x.len; j++) {
            if (j != i) {
                product *= (value - x.data[j]) / (x.data[i] - x.data[j]);
            }
        }
        sum += product;
    }
    return sum;
}

FUNC_DEF data_type linear_regression(const Vector x, const Vector y, data_type *const a, data_type *const b) {
    if ((x.len != y.len) || (a == NULL) || (b == NULL)) {
        return 0.0;  // Invalid operation
    }
    data_type sum_x = 0.0, sum_x_squared = 0.0;
    data_type sum_y = 0.0, sum_y_squared = 0.0;
    data_type sum_xy = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        sum_x += x.data[i];
        sum_x_squared += x.data[i] * x.data[i];
        sum_y += y.data[i];
        sum_y_squared += y.data[i] * y.data[i];
        sum_xy += x.data[i] * y.data[i];
    }
    *a = (x.len * sum_xy - sum_x * sum_y) / (x.len * sum_x_squared - sum_x * sum_x);
    *b = (sum_y * sum_x_squared - sum_xy * sum_x) / (x.len * sum_x_squared - sum_x * sum_x);
    data_type r = (x.len * sum_xy - sum_x * sum_y) / (square_root(x.len * sum_x_squared - sum_x * sum_x) * square_root(x.len * sum_y_squared - sum_y * sum_y));
    return r;
}

// Remember to free the returned vector after calling this function!
FUNC_DEF Vector polynomial_regression(const Vector x, const Vector y, const size_t order) {
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
                matrix_set(A, i, j, matrix_get(A, i, j) + power(x.data[l], (i + j)));
            }
            if (i != j) {
                matrix_set(A, j, i, matrix_get(A, i, j));
            }
        }
        for (size_t l = 0; l < x.len; l++) {
            b.data[i] += y.data[l] * power(x.data[l], i);
        }
    }
    Vector coefficients = gaussian_elimination(A, b);
    matrix_dealloc(&A);
    vector_dealloc(&b);
    return coefficients;
}

FUNC_DEF data_type compute_polynomial(const Vector coefficients, const data_type x) {
    data_type y = 0.0;
    for (size_t i = 0; i < coefficients.len; i++) {
        y += coefficients.data[i] * power(x, i);
    }
    return y;
}

#endif  // MATH_VECTOR_IMPLEMENTATION

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