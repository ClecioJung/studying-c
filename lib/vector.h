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

#ifndef __VECTOR_H
#define __VECTOR_H

// #include <math.h>
#include <stdbool.h>
// #include <stdint.h>
// #include <stdio.h>
#include <stdlib.h>
// #include <time.h>

typedef struct {
    double *data;
    size_t len;
} Vector;

Vector vector_alloc(const size_t len);
void vector_dealloc(Vector *const vector);
Vector vector_realloc(Vector *const vector, const size_t len);
bool vector_is_valid(const Vector vector);
void vector_set(const Vector vector, const size_t index, const double value);
double vector_get(const Vector vector, const size_t index);
Vector vector_random(const size_t len, const double min, const double max);
Vector vector_init(const size_t len, const double value);
void vector_assign(Vector *const vector, const Vector equals);
void vector_replace(Vector *const vector, const Vector equals);
Vector vector_copy(const Vector vector);
void vector_print(const Vector vector);
Vector vector_scale(const double scalar, const Vector vector);
double dot_product(const Vector a, const Vector b);
Vector cross_product(const Vector a, const Vector b);
double euclidean_norm(const Vector x);
double vector_max(const Vector x);
Vector vector_sum(const Vector a, const Vector b);
Vector vector_sub(const Vector a, const Vector b);
double vector_max_diff(Vector x, Vector previous_x);
bool vector_are_equal(const Vector a, const Vector b);

// 'over' methods override the contents of their arguments,
// avoiding the need to allocate more memory for the results
void vector_scale_over(const double scalar, const Vector *const vector);
void vector_copy_over(const Vector *const vector, const Vector vec_to_copy);

#endif  // __VECTOR_H

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