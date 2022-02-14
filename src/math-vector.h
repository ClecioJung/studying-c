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

#define MAX_ITERATIONS 10000
#define PRECISION 1e-10

typedef double data_type;

typedef struct {
    data_type *data;
    size_t len;
} Vector;

typedef struct {
    data_type **data;
    size_t rows;
    size_t cols;
} Matrix;

// Scalar functions
void swap(data_type *const a, data_type *const b);
uint64_t factorial(uint64_t value);
bool are_close(const data_type a, const data_type b, const data_type delta);
data_type maximum(const data_type a, const data_type b);
data_type minimum(const data_type a, const data_type b);
data_type generate_random_data_type(const data_type min, const data_type max);
data_type square_root(const data_type value);
data_type power(const data_type base, uint64_t expoent);
data_type exponential(const data_type value);

// Vector functions
Vector alloc_vector(const size_t len);
Vector realloc_vector(Vector *const vector, const size_t len);
void free_vector(Vector *const vector);
Vector random_vector(const size_t len, const data_type min, const data_type max);
Vector copy_vector(const Vector vector);
void print_vector(const Vector vector);
Vector scale_vector(const data_type scalar, const Vector vector);
data_type dot_product(const Vector a, const Vector b);
Vector cross_product(const Vector a, const Vector b);
data_type euclidean_norm(const Vector x);
data_type vector_max(const Vector x);
Vector sum_vectors(const Vector a, const Vector b);
Vector sub_vectors(const Vector a, const Vector b);
bool vectors_are_equal(const Vector a, const Vector b);

// Sorting functions
bool vector_is_sorted(const Vector vec);
void bubble_sort(const Vector vec);
void select_sort(const Vector vec);
void insert_sort(const Vector vec);
void shell_sort(const Vector vec);
void heap_sort(const Vector vec);
void merge_sort(const Vector vec);
void quicksort(const Vector vec);

// Search functions
size_t sequential_search(const Vector vec, const data_type value);
size_t binary_search(const Vector vec, const data_type value);

// Matrix functions
Matrix alloc_matrix(const size_t rows, const size_t cols);
void free_matrix(Matrix *const matrix);
Matrix random_matrix(const size_t rows, const size_t cols, const data_type min, const data_type max);
Matrix init_matrix(const size_t rows, const size_t cols, const data_type value);
Matrix identity_matrix(const size_t rows);
Matrix copy_matrix(const Matrix matrix);
void print_matrix(const Matrix matrix);
Matrix scale_matrix(const data_type scalar, const Matrix A);
Matrix mul_matrices(const Matrix A, const Matrix B);
Matrix mul_3_matrices(const Matrix A, const Matrix B, const Matrix C);
Vector mul_matrix_vector(const Matrix A, const Vector b);
Matrix matrix_transpose(const Matrix A);
Matrix matrix_symmetric(const Matrix A);
Matrix matrix_skew_symmetric(const Matrix A);
data_type trace(const Matrix A);
data_type determinant(const Matrix A);
void lu_decomposition(const Matrix A, Matrix *const L, Matrix *const U);
void lu_crout_decomposition(const Matrix A, Matrix *const L, Matrix *const U);
Matrix matrix_inverse(const Matrix A);
Matrix pseudo_inverse(const Matrix A);
Matrix sum_matrices(const Matrix A, const Matrix B);
Matrix sub_matrices(const Matrix A, const Matrix B);
bool matrices_are_equal(const Matrix A, const Matrix B);

// Methods for solving linear systems
Vector back_substitution(const Matrix A, const Vector b);
Vector forward_substitution(const Matrix A, const Vector b);
Vector gaussian_elimination(const Matrix A, const Vector b);
Vector gauss_jordan(const Matrix A, const Vector b);
Vector lu_solving(const Matrix A, const Vector b);
data_type get_max_diff(Vector x, Vector previous_x);
Vector jacobi_method(const Matrix A, const Vector b);
Vector gauss_seidel(const Matrix A, const Vector b);
bool columns_condition(const Matrix A);
bool rows_condition(const Matrix A);
bool sassenfeld_condition(const Matrix A);

// Curve fitting functions
data_type lagrange_interpolation(const Vector x, const Vector y, const data_type value);
data_type linear_regression(const Vector x, const Vector y, data_type *const a, data_type *const b);
Vector polynomial_regression(const Vector x, const Vector y, const size_t order);
data_type compute_polynomial(const Vector coefficients, const data_type x);

#endif  // __MATH_VECTOR_H

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------

#ifdef MATH_VECTOR_IMPLEMENTATION

void swap(data_type *const a, data_type *const b) {
    if ((a != NULL) && (b != NULL)) {
        const data_type temp = *a;
        *a = *b;
        *b = temp;
    }
}

uint64_t factorial(uint64_t value) {
    uint64_t result = 1;
    while (value > 1) {
        result = result * value;
        value--;
    }
    return result;
}

bool are_close(const data_type a, const data_type b, const data_type delta) {
    return (fabs(a - b) < delta);
}

data_type maximum(const data_type a, const data_type b) {
    return ((a > b) ? a : b);
}

data_type minimum(const data_type a, const data_type b) {
    return ((a < b) ? a : b);
}

data_type generate_random_data_type(const data_type min, const data_type max) {
    const data_type rand_unitary = ((data_type)rand()) / ((data_type)RAND_MAX);
    return rand_unitary * (max - min) + min;
}

// My own square root function, so I don't need to link with -lm,
// avoiding any dependencies
data_type square_root(const data_type value) {
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
data_type power(const data_type base, uint64_t expoent) {
    data_type result = 1.0;
    while (expoent != 0) {
        result *= base;
        expoent--;
    }
    return result;
}

// My own exponential function, so I don't need to link with -lm,
// avoiding any dependencies
data_type exponential(const data_type value) {
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

// Remember to free the vector after calling this function!
Vector alloc_vector(const size_t len) {
    Vector vector = (Vector){
        .data = (data_type *)malloc(len * sizeof(data_type)),
        .len = len,
    };
    if (vector.data == NULL) {
        vector.len = 0;
    }
    return vector;
}

// Remember to free the vector after calling this function!
Vector realloc_vector(Vector *const vector, const size_t len) {
    data_type *new_data = (data_type *)realloc(vector->data, len * sizeof(data_type));
    if (new_data != NULL) {
        vector->len = len;
        vector->data = new_data;
    }
    return *vector;
}

void free_vector(Vector *const vector) {
    if (vector == NULL) {
        return;
    }
    free(vector->data);
    vector->data = NULL;
    vector->len = 0;
}

// Remember to free the vector after calling this function!
Vector random_vector(const size_t len, const data_type min, const data_type max) {
    Vector vector = alloc_vector(len);
    if (vector.data != NULL) {
        srand(time(NULL));
        for (size_t i = 0; i < len; i++) {
            vector.data[i] = generate_random_data_type(min, max);
        }
    }
    return vector;
}

// Remember to free the vector after calling this function!
Vector init_vector(const size_t len, const data_type value) {
    Vector vector = alloc_vector(len);
    if (vector.data != NULL) {
        for (size_t i = 0; i < len; i++) {
            vector.data[i] = value;
        }
    }
    return vector;
}

// Remember to free the vector after calling this function!
Vector copy_vector(const Vector vector) {
    Vector new_vec = alloc_vector(vector.len);
    if (new_vec.data != NULL) {
        for (size_t i = 0; i < vector.len; i++) {
            new_vec.data[i] = vector.data[i];
        }
    }
    return new_vec;
}

void print_vector(const Vector vector) {
    for (size_t i = 0; i < vector.len; i++) {
        printf("[%03ld]: %g\n", i, vector.data[i]);
    }
    printf("\n");
}

// Remember to free the vector after calling this function!
Vector scale_vector(const data_type scalar, const Vector vector) {
    Vector new_vec = alloc_vector(vector.len);
    if (new_vec.data != NULL) {
        for (size_t i = 0; i < vector.len; i++) {
            new_vec.data[i] = scalar * vector.data[i];
        }
    }
    return new_vec;
}

data_type dot_product(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return 0.0;  // Invalid operation
    }
    data_type result = 0.0;
    for (size_t i = 0; i < a.len; i++) {
        result += a.data[i] * b.data[i];
    }
    return result;
}

// Remember to free the vector after calling this function!
Vector cross_product(const Vector a, const Vector b) {
    if ((a.len != 3) || (b.len != 3)) {
        return (Vector){0};  // Invalid operation
    }
    Vector result = alloc_vector(a.len);
    result.data[0] = a.data[1] * b.data[2] - a.data[2] * b.data[1];
    result.data[1] = a.data[2] * b.data[0] - a.data[0] * b.data[2];
    result.data[2] = a.data[0] * b.data[1] - a.data[1] * b.data[0];
    return result;
}

data_type euclidean_norm(const Vector x) {
    data_type value = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        value += x.data[i] * x.data[i];
    }
    return square_root(value);
}

data_type vector_max(const Vector x) {
    data_type value = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        value = maximum(value, x.data[i]);
    }
    return value;
}

// Remember to free the vector after calling this function!
Vector sum_vectors(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector new_vec = alloc_vector(a.len);
    if (new_vec.data != NULL) {
        for (size_t i = 0; i < a.len; i++) {
            new_vec.data[i] = a.data[i] + b.data[i];
        }
    }
    return new_vec;
}

// Remember to free the vector after calling this function!
Vector sub_vectors(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector new_vec = alloc_vector(a.len);
    if (new_vec.data != NULL) {
        for (size_t i = 0; i < a.len; i++) {
            new_vec.data[i] = a.data[i] - b.data[i];
        }
    }
    return new_vec;
}

bool vectors_are_equal(const Vector a, const Vector b) {
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

bool vector_is_sorted(const Vector vec) {
    for (size_t i = 0; (i + 1) < vec.len; i++) {
        if (vec.data[i] > vec.data[i + 1]) {
            return false;
        }
    }
    return true;
}

void bubble_sort(const Vector vec) {
    for (size_t i = 1; i < vec.len; i++) {
        for (size_t j = (vec.len - 1); j >= i; j--) {
            if (vec.data[j - 1] > vec.data[j]) {
                swap(&vec.data[j - 1], &vec.data[j]);
            }
        }
    }
}

void select_sort(const Vector vec) {
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

void insert_sort(const Vector vec) {
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

void shell_sort(const Vector vec) {
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
void heapify(data_type *const array, const size_t len, const size_t index) {
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
void heap_sort(const Vector vec) {
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
void merge(data_type *const array, const size_t p, const size_t q, const size_t r) {
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
void split_and_sort(data_type *const array, const size_t left, const size_t right) {
    if (left < right) {
        data_type middle = left + (right - left) / 2;
        split_and_sort(array, left, middle);
        split_and_sort(array, middle + 1, right);
        merge(array, left, middle, right);
    }
}
void merge_sort(const Vector vec) {
    if (vec.len > 0) {
        split_and_sort(vec.data, 0, vec.len - 1);
    }
}

void qs(data_type *const array, const size_t left, const size_t right) {
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
void quicksort(const Vector vec) {
    if (vec.len > 0) {
        qs(vec.data, 0, vec.len - 1);
    }
}

size_t sequential_search(const Vector vec, const data_type value) {
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
size_t binary_search(const Vector vec, const data_type value) {
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
Matrix alloc_matrix(const size_t rows, const size_t cols) {
    Matrix matrix = (Matrix){
        .data = (data_type **)malloc(rows * sizeof(data_type *)),
        .rows = rows,
        .cols = cols,
    };
    if (matrix.data != NULL) {
        for (size_t i = 0; i < rows; i++) {
            matrix.data[i] = NULL;
        }
        for (size_t i = 0; i < rows; i++) {
            matrix.data[i] = (data_type *)malloc(cols * sizeof(data_type));
            if (matrix.data[i] == NULL) {
                free_matrix(&matrix);
                break;
            }
        }
    } else {
        matrix.rows = 0;
        matrix.cols = 0;
    }
    return matrix;
}

void free_matrix(Matrix *const matrix) {
    if (matrix == NULL) {
        return;
    }
    for (size_t i = 0; i < matrix->rows; i++) {
        free(matrix->data[i]);
    }
    free(matrix->data);
    matrix->data = NULL;
    matrix->rows = 0;
    matrix->cols = 0;
}

// Remember to free the matrix after calling this function!
Matrix random_matrix(const size_t rows, const size_t cols, const data_type min, const data_type max) {
    Matrix matrix = alloc_matrix(rows, cols);
    if (matrix.data != NULL) {
        srand(time(NULL));
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                matrix.data[i][j] = generate_random_data_type(min, max);
            }
        }
    }
    return matrix;
}

// Remember to free the matrix after calling this function!
Matrix init_matrix(const size_t rows, const size_t cols, const data_type value) {
    Matrix matrix = alloc_matrix(rows, cols);
    if (matrix.data != NULL) {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < cols; j++) {
                matrix.data[i][j] = value;
            }
        }
    }
    return matrix;
}

// Remember to free the matrix after calling this function!
Matrix identity_matrix(const size_t rows) {
    Matrix matrix = alloc_matrix(rows, rows);
    if (matrix.data != NULL) {
        for (size_t i = 0; i < rows; i++) {
            for (size_t j = 0; j < rows; j++) {
                matrix.data[i][j] = (i == j ? 1 : 0);
            }
        }
    }
    return matrix;
}

// Remember to free the matrix after calling this function!
Matrix copy_matrix(const Matrix matrix) {
    Matrix new_matrix = alloc_matrix(matrix.rows, matrix.cols);
    if (new_matrix.data != NULL) {
        for (size_t i = 0; i < matrix.rows; i++) {
            for (size_t j = 0; j < matrix.cols; j++) {
                new_matrix.data[i][j] = matrix.data[i][j];
            }
        }
    }
    return new_matrix;
}

void print_matrix(const Matrix matrix) {
    printf("%*s", 2, "");
    for (size_t j = 0; j < matrix.cols; j++) {
        printf("%*s[%03ld] ", 5, "", j);
    }
    printf("\n");
    for (size_t i = 0; i < matrix.rows; i++) {
        printf("[%03ld]: ", i);
        for (size_t j = 0; j < matrix.cols; j++) {
            printf("%-10g ", matrix.data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Remember to free the vector after calling this function!
Matrix scale_matrix(const data_type scalar, const Matrix A) {
    Matrix result = alloc_matrix(A.rows, A.cols);
    if (result.data != NULL) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                result.data[i][j] = scalar * A.data[i][j];
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix mul_matrices(const Matrix A, const Matrix B) {
    if (A.cols != B.rows) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = init_matrix(A.rows, B.cols, 0.0);
    if (result.data != NULL) {
        for (size_t i = 0; i < result.rows; i++) {
            for (size_t j = 0; j < result.cols; j++) {
                for (size_t k = 0; k < A.cols; k++) {
                    result.data[i][j] += A.data[i][k] * B.data[k][j];
                }
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix mul_3_matrices(const Matrix A, const Matrix B, const Matrix C) {
    if ((A.cols != B.rows) || (B.cols != C.rows)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix mul = mul_matrices(A, B);
    Matrix result = mul_matrices(mul, C);
    free_matrix(&mul);
    return result;
}

// Remember to free the returned vector after calling this function!
Vector mul_matrix_vector(const Matrix A, const Vector b) {
    if (A.cols != b.len) {
        return (Vector){0};  // Invalid operation
    }
    Vector result = init_vector(A.rows, 0.0);
    if (result.data != NULL) {
        for (size_t i = 0; i < result.len; i++) {
            for (size_t k = 0; k < A.cols; k++) {
                result.data[i] += A.data[i][k] * b.data[k];
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_transpose(const Matrix A) {
    Matrix transpose = alloc_matrix(A.cols, A.rows);
    if (transpose.data != NULL) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                transpose.data[j][i] = A.data[i][j];
            }
        }
    }
    return transpose;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_symmetric(const Matrix A) {
    if (A.rows != A.cols) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix sym = alloc_matrix(A.rows, A.cols);
    if (sym.data != NULL) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                sym.data[i][j] = (A.data[i][j] + A.data[j][i])/2.0;
            }
        }
    }
    return sym;
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_skew_symmetric(const Matrix A) {
    if (A.rows != A.cols) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix skew = alloc_matrix(A.rows, A.cols);
    if (skew.data != NULL) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                skew.data[i][j] = (A.data[i][j] - A.data[j][i])/2.0;
            }
        }
    }
    return skew;
}

data_type trace(const Matrix A) {
    if (A.rows != A.cols) {
        return 0.0;  // Invalid operation
    }
    data_type trace = 0.0;
    for (size_t i = 0; i < A.rows; i++) {
        trace += A.data[i][i];
    }
    return trace;
}

data_type determinant(const Matrix A) {
    if (A.rows != A.cols) {
        return 0.0;  // Invalid operation
    }
    // Gaussian elimination
    Matrix B = copy_matrix(A);
    if (B.data == NULL) {
        return 0.0;
    }
    for (size_t k = 0; (k + 1) < B.rows; k++) {
        for (size_t i = (k + 1); i < B.rows; i++) {
            const data_type m = B.data[i][k] / B.data[k][k];
            for (size_t j = (k + 1); j < B.rows; j++) {
                B.data[i][j] -= m * B.data[k][j];
            }
        }
    }
    data_type det = 1.0;
    for (size_t k = 0; k < B.rows; k++) {
        det *= B.data[k][k];
    }
    free_matrix(&B);
    return det;
}

// Remember to free L and U matrices after calling this function!
void lu_decomposition(const Matrix A, Matrix *const L, Matrix *const U) {
    if ((A.rows != A.cols) || L == NULL || U == NULL) {
        return;  // Invalid operation
    }
    *L = identity_matrix(A.rows);
    *U = copy_matrix(A);
    if ((L->data == NULL) || (U->data == NULL)) {
        free_matrix(L);
        free_matrix(U);
        return;
    }
    for (size_t k = 0; (k + 1) < U->rows; k++) {
        for (size_t i = (k + 1); i < U->rows; i++) {
            L->data[i][k] = U->data[i][k] / U->data[k][k];
            U->data[i][k] = 0.0;
            for (size_t j = (k + 1); j < U->rows; j++) {
                U->data[i][j] -= L->data[i][k] * U->data[k][j];
            }
        }
    }
}

// Remember to free L and U matrices after calling this function!
void lu_crout_decomposition(const Matrix A, Matrix *const L, Matrix *const U) {
    if ((A.rows != A.cols) || L == NULL || U == NULL) {
        return;  // Invalid operation
    }
    // Crout Decomposition
    *L = init_matrix(A.rows, A.cols, 0.0);
    *U = identity_matrix(A.rows);
    if ((L->data == NULL) || (U->data == NULL)) {
        free_matrix(L);
        free_matrix(U);
        return;
    }
    for (size_t k = 0; k < A.rows; k++) {
        for (size_t i = k; i < A.rows; i++) {
            L->data[i][k] = A.data[i][k];
            for (size_t l = 0; l < k; l++) {
                L->data[i][k] -= L->data[i][l] * U->data[l][k];
            }
        }
        if ((k + 1) != A.rows) {
            for (size_t j = (k + 1); j < A.rows; j++) {
                U->data[k][j] = A.data[k][j];
                for (size_t l = 0; l < k; l++) {
                    U->data[k][j] -= L->data[k][l] * U->data[l][j];
                }
                U->data[k][j] /= L->data[k][k];
            }
        }
    }
}

// Remember to free the returned matrix after calling this function!
Matrix matrix_inverse(const Matrix A) {
    if (A.rows != A.cols) {
        return (Matrix){0};  // Invalid operation
    }
    // Computes the inverse matrix of A using LU decomposition
    Matrix L, U;
    lu_decomposition(A, &L, &U);
    Matrix inv = alloc_matrix(A.rows, A.cols);
    Vector b = alloc_vector(A.rows);
    if ((L.data == NULL) || (U.data == NULL) || (inv.data == NULL) || (b.data == NULL)) {
        free_matrix(&L);
        free_matrix(&U);
        free_matrix(&inv);
        free_vector(&b);
        return (Matrix){0};
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.rows; j++) {
            b.data[j] = (j == i) ? 1.0 : 0.0;
        }
        Vector d = forward_substitution(L, b);
        if (d.data == NULL) {
            free_matrix(&inv);
            break;
        }
        Vector x = back_substitution(U, d);
        if (x.data == NULL) {
            free_vector(&d);
            free_matrix(&inv);
            break;
        }
        for (size_t j = 0; j < A.rows; j++) {
            inv.data[j][i] = x.data[j];
        }
        free_vector(&d);
        free_vector(&x);
    }
    free_matrix(&L);
    free_matrix(&U);
    free_vector(&b);
    return inv;
}

Matrix pseudo_inverse(const Matrix A) {
    Matrix pseudo_inv, mul, inv;
    Matrix transpose = matrix_transpose(A);
    if (A.rows > A.cols) {
        // Left pseudo inverse
        // A^+ = (A^T * A)^-1 * A^T
        mul = mul_matrices(transpose, A);
        inv = matrix_inverse(mul);
        pseudo_inv = mul_matrices(inv, transpose);
    } else {
        // Right pseudo inverse
        // A^+ = A^T * (A * A^T)^-1
        mul = mul_matrices(A, transpose);
        inv = matrix_inverse(mul);
        pseudo_inv = mul_matrices(transpose, inv);
    }
    free_matrix(&mul);
    free_matrix(&inv);
    free_matrix(&transpose);
    return pseudo_inv;
}

// Remember to free the returned matrix after calling this function!
Matrix sum_matrices(const Matrix A, const Matrix B) {
    if ((A.rows != B.rows) || (A.cols != B.cols)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = alloc_matrix(A.rows, A.cols);
    if (result.data != NULL) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                result.data[i][j] = A.data[i][j] + B.data[i][j];
            }
        }
    }
    return result;
}

// Remember to free the returned matrix after calling this function!
Matrix sub_matrices(const Matrix A, const Matrix B) {
    if ((A.rows != B.rows) || (A.cols != B.cols)) {
        return (Matrix){0};  // Invalid operation
    }
    Matrix result = alloc_matrix(A.rows, A.cols);
    if (result.data != NULL) {
        for (size_t i = 0; i < A.rows; i++) {
            for (size_t j = 0; j < A.cols; j++) {
                result.data[i][j] = A.data[i][j] - B.data[i][j];
            }
        }
    }
    return result;
}

bool matrices_are_equal(const Matrix A, const Matrix B) {
    if ((A.rows != B.rows) || (A.cols != B.cols)) {
        return false;
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            if (!are_close(A.data[i][j], B.data[i][j], PRECISION)) {
                return false;
            }
        }
    }
    return true;
}

// Remember to free the returned vector after calling this function!
Vector back_substitution(const Matrix A, const Vector b) {
    if ((A.rows != A.cols) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = alloc_vector(b.len);
    if (x.data != NULL) {
        for (size_t i = (b.len - 1); i < b.len; i--) {
            data_type sum = b.data[i];
            for (size_t j = (i + 1); j < b.len; j++) {
                sum -= A.data[i][j] * x.data[j];
            }
            x.data[i] = sum / A.data[i][i];
        }
    }
    return x;
}

// Remember to free the returned vector after calling this function!
Vector forward_substitution(const Matrix A, const Vector b) {
    if ((A.rows != A.cols) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = alloc_vector(b.len);
    if (x.data != NULL) {
        for (size_t i = 0; i < b.len; i++) {
            data_type sum = b.data[i];
            for (size_t j = 0; j < i; j++) {
                sum -= A.data[i][j] * x.data[j];
            }
            x.data[i] = sum / A.data[i][i];
        }
    }
    return x;
}

// Remember to free the returned vector after calling this function!
Vector gaussian_elimination(const Matrix A, const Vector b) {
    if ((A.rows != A.cols) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Matrix A_copied = copy_matrix(A);
    Vector b_copied = copy_vector(b);
    if ((A_copied.data == NULL) || (b_copied.data == NULL)) {
        return (Vector){0};
    }
    // Partial pivoting
    for (size_t k = 0; (k + 1) < A_copied.rows; k++) {
        data_type w = fabs(A_copied.data[k][k]);
        size_t r = k;
        for (size_t i = (k + 1); i < A_copied.rows; i++) {
            if (fabs(A_copied.data[i][k]) > w) {
                w = fabs(A_copied.data[i][k]);
                r = i;
            }
        }
        if (r != k) {
            for (size_t i = k; i < A_copied.rows; i++) {
                w = A_copied.data[k][i];
                A_copied.data[k][i] = A_copied.data[r][i];
                A_copied.data[r][i] = w;
            }
            w = b_copied.data[k];
            b_copied.data[k] = b_copied.data[r];
            b_copied.data[r] = w;
        }
    }
    // Gaussian elimination
    for (size_t k = 0; (k + 1) < A_copied.rows; k++) {
        for (size_t i = (k + 1); i < A_copied.rows; i++) {
            const data_type m = A_copied.data[i][k] / A_copied.data[k][k];
            for (size_t j = (k + 1); j < A_copied.rows; j++) {
                A_copied.data[i][j] -= m * A_copied.data[k][j];
            }
            b_copied.data[i] -= m * b_copied.data[k];
        }
    }
    Vector x = back_substitution(A_copied, b_copied);
    free_matrix(&A_copied);
    free_vector(&b_copied);
    return x;
}

// Remember to free the returned vector after calling this function!
Vector gauss_jordan(const Matrix A, const Vector b) {
    if ((A.rows != A.cols) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Matrix A_copied = copy_matrix(A);
    Vector b_copied = copy_vector(b);
    if ((A_copied.data == NULL) || (b_copied.data == NULL)) {
        return (Vector){0};
    }
    for (size_t k = 0; k < A_copied.rows; k++) {
        for (size_t i = 0; i < A_copied.rows; i++) {
            data_type m = A_copied.data[i][k] / A_copied.data[k][k];
            data_type p = A_copied.data[k][k];
            if (i == k) {
                for (size_t j = k; j < A_copied.cols; j++) {
                    A_copied.data[i][j] /= p;
                }
                b_copied.data[i] /= p;
            } else {
                for (size_t j = k; j < A_copied.cols; j++) {
                    A_copied.data[i][j] -= m * A_copied.data[k][j];
                }
                b_copied.data[i] -= m * b_copied.data[k];
            }
        }
    }
    free_matrix(&A_copied);
    return b_copied;
}

// Remember to free the returned vector after calling this function!
Vector lu_solving(const Matrix A, const Vector b) {
    if ((A.rows != A.cols) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = (Vector){0};
    Matrix L, U;
    lu_decomposition(A, &L, &U);
    if ((L.data != NULL) && (U.data != NULL)) {
        Vector d = forward_substitution(L, b);
        if (d.data != NULL) {
            x = back_substitution(U, d);
        }
        free_vector(&d);
    }
    free_matrix(&L);
    free_matrix(&U);
    return x;
}

data_type get_max_diff(const Vector x, const Vector previous_x) {
    if (x.len != previous_x.len) {
        return INFINITY;  // Invalid operation
    }
    data_type error = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        error = maximum(fabs(x.data[i] - previous_x.data[i]), error);
    }
    return error;
}

// Remember to free the returned vector after calling this function!
Vector jacobi_method(const Matrix A, const Vector b) {
    if ((A.rows != A.cols) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = init_vector(b.len, 0.0);
    if (x.data != NULL) {
        for (size_t k = 0; k < MAX_ITERATIONS; k++) {
            Vector previous_x = copy_vector(x);
            if (previous_x.data == NULL) {
                free_vector(&x);
                break;
            }
            for (size_t i = 0; i < x.len; i++) {
                data_type sum = b.data[i];
                for (size_t j = 0; j < x.len; j++) {
                    if (i != j) {
                        sum -= A.data[i][j] * previous_x.data[j];
                    }
                }
                x.data[i] = sum / A.data[i][i];
            }
            const data_type error = get_max_diff(x, previous_x);
            free_vector(&previous_x);
            if (error < PRECISION) {
                break;
            }
        }
    }
    return x;
}

// Remember to free the returned vector after calling this function!
Vector gauss_seidel(const Matrix A, const Vector b) {
    if ((A.rows != A.cols) || (A.rows != b.len)) {
        return (Vector){0};  // Invalid operation
    }
    Vector x = init_vector(b.len, 0.0);
    if (x.data != NULL) {
        for (size_t k = 0; k < MAX_ITERATIONS; k++) {
            Vector previous_x = copy_vector(x);
            if (previous_x.data == NULL) {
                free_vector(&x);
                break;
            }
            for (size_t i = 0; i < b.len; i++) {
                data_type sum = b.data[i];
                for (size_t j = 0; j < x.len; j++) {
                    if (i != j) {
                        sum -= A.data[i][j] * x.data[j];
                    }
                }
                x.data[i] = sum / A.data[i][i];
            }
            const data_type error = get_max_diff(x, previous_x);
            free_vector(&previous_x);
            if (error < PRECISION) {
                break;
            }
        }
    }
    return x;
}

bool columns_condition(const Matrix A) {
    if (A.rows != A.cols) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        data_type sum = 0.0;
        for (size_t j = 0; j < A.rows; j++) {
            if (i != j) {
                sum += fabs(A.data[j][i]);
            }
        }
        if (fabs(A.data[i][i]) < sum) {
            return false;
        }
    }
    return true;
}

bool rows_condition(const Matrix A) {
    if (A.rows != A.cols) {
        return false;  // Invalid operation
    }
    for (size_t i = 0; i < A.rows; i++) {
        data_type sum = 0.0;
        for (size_t j = 0; j < A.rows; j++) {
            if (i != j) {
                sum += fabs(A.data[i][j]);
            }
        }
        if (fabs(A.data[i][i]) < sum) {
            return false;
        }
    }
    return true;
}

bool sassenfeld_condition(const Matrix A) {
    if (A.rows != A.cols) {
        return false;  // Invalid operation
    }
    bool condition = false;
    Vector beta = alloc_vector(A.rows);
    if (beta.data != NULL) {
        for (size_t i = 0; i < A.rows; i++) {
            data_type sum = 0.0;
            for (size_t j = 0; j < i; j++) {
                sum += fabs(A.data[i][j]) * beta.data[j];
            }
            for (size_t j = (i + 1); j < A.rows; j++) {
                sum += fabs(A.data[i][j]);
            }
            beta.data[i] = sum / A.data[i][i];
        }
        if (vector_max(beta) < 1.0) {
            condition = true;
        }
        free_vector(&beta);
    }
    return condition;
}

data_type lagrange_interpolation(const Vector x, const Vector y, const data_type value) {
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

data_type linear_regression(const Vector x, const Vector y, data_type *const a, data_type *const b) {
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
Vector polynomial_regression(const Vector x, const Vector y, const size_t order) {
    if ((x.len != y.len) || (x.len < (order + 1))) {
        return (Vector){0};  // Invalid operation
    }
    Matrix A = init_matrix(order + 1, order + 1, 0.0);
    Vector b = init_vector(order + 1, 0.0);
    if ((A.data == NULL) || (b.data == NULL)) {
        free_matrix(&A);
        free_vector(&b);
        return (Vector){0};
    }
    // Computing augmented matrix
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j <= i; j++) {
            for (size_t l = 0; l < x.len; l++) {
                A.data[i][j] += power(x.data[l], (i + j));
            }
            if (i != j) {
                A.data[j][i] = A.data[i][j];
            }
        }
        for (size_t l = 0; l < x.len; l++) {
            b.data[i] += y.data[l] * power(x.data[l], i);
        }
    }
    Vector coefficients = gaussian_elimination(A, b);
    free_matrix(&A);
    free_vector(&b);
    return coefficients;
}

data_type compute_polynomial(const Vector coefficients, const data_type x) {
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