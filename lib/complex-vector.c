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

#include "complex-vector.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>

#include "complex.h"
#include "scalar.h"
#include "vector.h"

#define MAX_ITERATIONS 10000
#define COMPARATION_PRECISION 1e-8

// Remember to free the returned vector after calling this function!
Complex_Vector complex_vector_alloc(const size_t len) {
    Complex_Vector vector = (Complex_Vector){0};
    if (len != 0) {
        vector.data = (Complex *)malloc(len * sizeof(Complex));
        if (vector.data != NULL) {
            vector.len = len;
        }
    }
    return vector;
}

void complex_vector_dealloc(Complex_Vector *const vector) {
    if (vector == NULL) {
        return;
    }
    free(vector->data);
    vector->data = NULL;
    vector->len = 0;
}

bool complex_vector_is_valid(const Complex_Vector vector) {
    return (vector.data != NULL);
}

void complex_vector_set(const Complex_Vector vector, const size_t index, const Complex value) {
    if (!complex_vector_is_valid(vector) || (index >= vector.len)) {
        return;
    }
    vector.data[index] = value;
}

Complex complex_vector_get(const Complex_Vector vector, const size_t index) {
    if (!complex_vector_is_valid(vector) || (index >= vector.len)) {
        return complex_init(NAN, NAN);
    }
    return vector.data[index];
}

// Remember to free the returned vector after calling this function!
Complex_Vector complex_vector_init(const size_t len, const Complex value) {
    Complex_Vector vector = complex_vector_alloc(len);
    if (complex_vector_is_valid(vector)) {
        complex_vector_init_over(vector, value);
    }
    return vector;
}

void complex_vector_print(const Complex_Vector vector) {
    if (complex_vector_is_valid(vector)) {
        for (size_t i = 0; i < vector.len; i++) {
            printf("[%03zu]: ", i);
            complex_print(vector.data[i]);
            printf("\n");
        }
        printf("\n");
    }
}

// Remember to free the returned vector after calling this function!
Vector complex_vector_get_real_part(const Complex_Vector vector) {
    Vector real_vector = vector_alloc(vector.len);
    if (vector_is_valid(real_vector)) {
        for (size_t i = 0; i < real_vector.len; i++) {
            real_vector.data[i] = vector.data[i].real;
        }
    }
    return real_vector;
}

// Remember to free the returned vector after calling this function!
Vector complex_vector_get_imag_part(const Complex_Vector vector) {
    Vector imag_vector = vector_alloc(vector.len);
    if (vector_is_valid(imag_vector)) {
        for (size_t i = 0; i < imag_vector.len; i++) {
            imag_vector.data[i] = vector.data[i].imag;
        }
    }
    return imag_vector;
}

void complex_vector_init_over(const Complex_Vector vector, const Complex value) {
    if (complex_vector_is_valid(vector)) {
        for (size_t i = 0; i < vector.len; i++) {
            vector.data[i] = value;
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