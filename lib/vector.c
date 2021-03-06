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

#include "vector.h"

#include <math.h>
#include <stdarg.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "scalar.h"

#define MAX_ITERATIONS 10000
#define COMPARATION_PRECISION 1e-8

// Remember to free the returned vector after calling this function!
Vector vector_alloc(const size_t len) {
    Vector vector = (Vector){0};
    if (len != 0) {
        vector.data = (double *)malloc(len * sizeof(double));
        if (vector.data != NULL) {
            vector.len = len;
        }
    }
    return vector;
}

void vector_dealloc(Vector *const vector) {
    if (vector == NULL) {
        return;
    }
    free(vector->data);
    vector->data = NULL;
    vector->len = 0;
}

// Remember to free the returned vector after calling this function!
Vector vector_realloc(Vector *const vector, const size_t len) {
    if (len != 0) {
        double *new_data = (double *)realloc(vector->data, len * sizeof(double));
        if (new_data != NULL) {
            vector->len = len;
            vector->data = new_data;
        }
    } else {
        vector_dealloc(vector);
    }
    return *vector;
}

bool vector_is_valid(const Vector vector) {
    return (vector.data != NULL);
}

void vector_set(const Vector vector, const size_t index, const double value) {
    if (!vector_is_valid(vector) || (index >= vector.len)) {
        return;
    }
    vector.data[index] = value;
}

double vector_get(const Vector vector, const size_t index) {
    if (!vector_is_valid(vector) || (index >= vector.len)) {
        return NAN;
    }
    return vector.data[index];
}

// Remember to free the returned vector after calling this function!
Vector vector_random(const size_t len, const double min, const double max) {
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
Vector vector_init(const size_t len, const double value) {
    Vector vector = vector_alloc(len);
    if (vector_is_valid(vector)) {
        vector_init_over(vector, value);
    }
    return vector;
}

// Remember to free the returned vector after calling this function!
Vector vector_copy(const Vector vector) {
    Vector new_vec = vector_alloc(vector.len);
    if (vector_is_valid(new_vec)) {
        vector_copy_over(new_vec, vector);
    }
    return new_vec;
}

// Remember to free the returned vector after calling this function!
// Incorrect calls to this function will not result in compilation errors.
// Variable arguments are unsafe by design.
Vector vector_new(const size_t len, ...) {
    Vector vector = vector_alloc(len);
    if (vector_is_valid(vector)) {
        va_list valist;
        va_start(valist, len);
        for (size_t i = 0; i < len; i++) {
            vector.data[i] = (double)va_arg(valist, double);
        }
        va_end(valist);
    }
    return vector;
}

void vector_print(const Vector vector) {
    for (size_t i = 0; i < vector.len; i++) {
        const double value = fabs(vector.data[i]) > COMPARATION_PRECISION ? vector.data[i] : 0.0;
        printf("[%03zu]: %lg\n", i, value);
    }
    printf("\n");
}

// Remember to free the returned vector after calling this function!
Vector vector_scale(const double scalar, const Vector vector) {
    Vector new_vec = vector_alloc(vector.len);
    if (vector_is_valid(new_vec)) {
        for (size_t i = 0; i < vector.len; i++) {
            new_vec.data[i] = scalar * vector.data[i];
        }
    }
    return new_vec;
}

double vector_dot_product(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return 0.0;  // Invalid operation
    }
    double result = 0.0;
    for (size_t i = 0; i < a.len; i++) {
        result += a.data[i] * b.data[i];
    }
    return result;
}

// Remember to free the returned vector after calling this function!
Vector vector_cross_product(const Vector a, const Vector b) {
    if ((a.len != 3) || (b.len != 3)) {
        return (Vector){0};  // Invalid operation
    }
    Vector result = vector_alloc(a.len);
    result.data[0] = a.data[1] * b.data[2] - a.data[2] * b.data[1];
    result.data[1] = a.data[2] * b.data[0] - a.data[0] * b.data[2];
    result.data[2] = a.data[0] * b.data[1] - a.data[1] * b.data[0];
    return result;
}

double vector_norm(const Vector x) {
    double value = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        value += x.data[i] * x.data[i];
    }
    return square_root(value);
}

double vector_max(const Vector x) {
    double value = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        value = maximum(value, x.data[i]);
    }
    return value;
}

double vector_min(const Vector x) {
    double value = INFINITY;
    for (size_t i = 0; i < x.len; i++) {
        value = minimum(value, x.data[i]);
    }
    return value;
}

double vector_max_abs(const Vector x) {
    double error = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        error = maximum(fabs(x.data[i]), error);
    }
    return error;
}

double vector_mean(const Vector x) {
    double sum = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        sum += x.data[i];
    }
    return (sum / x.len);
}

// Remember to free the returned vector after calling this function!
Vector vector_sum(const Vector a, const Vector b) {
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
Vector vector_sub(const Vector a, const Vector b) {
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

double vector_max_diff(const Vector x, const Vector previous_x) {
    if (x.len != previous_x.len) {
        return INFINITY;  // Invalid operation
    }
    double error = 0.0;
    for (size_t i = 0; i < x.len; i++) {
        error = maximum(fabs(x.data[i] - previous_x.data[i]), error);
    }
    return error;
}

bool vector_are_equal(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return false;
    }
    for (size_t i = 0; i < a.len; i++) {
        if (!are_close(a.data[i], b.data[i], COMPARATION_PRECISION)) {
            return false;
        }
    }
    return true;
}

bool vector_is_null(const Vector vec) {
    for (size_t i = 0; i < vec.len; i++) {
        if (!are_close(vec.data[i], 0.0, COMPARATION_PRECISION)) {
            return false;
        }
    }
    return true;
}

bool vector_are_orthogonal(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return false;
    }
    return are_close(vector_dot_product(a, b), 0.0, COMPARATION_PRECISION);
}

void vector_init_over(const Vector vector, const double value) {
    for (size_t i = 0; i < vector.len; i++) {
        vector.data[i] = value;
    }
}

void vector_copy_over(const Vector vector, const Vector to_copy) {
    if (!vector_is_valid(vector) || (vector.len != to_copy.len)) {
        return;
    }
    for (size_t i = 0; i < vector.len; i++) {
        vector.data[i] = to_copy.data[i];
    }
}

void vector_scale_over(const double scalar, const Vector vector) {
    for (size_t i = 0; i < vector.len; i++) {
        vector.data[i] *= scalar;
    }
}

void vector_sum_over(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return;  // Invalid operation
    }
    for (size_t i = 0; i < a.len; i++) {
        a.data[i] += b.data[i];
    }
}

void vector_sub_over(const Vector a, const Vector b) {
    if (a.len != b.len) {
        return;  // Invalid operation
    }
    for (size_t i = 0; i < a.len; i++) {
        a.data[i] -= b.data[i];
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