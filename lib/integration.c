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

#include "integration.h"

#include <math.h>
#include <stdlib.h>

#include "linear-systems.h"
#include "matrix.h"
#include "polynomial.h"
#include "scalar.h"
#include "vector.h"

double trapezoidal_rule(const integration_fn fn, const double start, const double end, const size_t segments) {
    double step = (end - start) / segments;
    double sum = (fn(start) + fn(end)) / 2.0;
    for (size_t i = 1; i < segments; i++) {
        sum += fn(start + i * step);
    }
    return (step * sum);
}

double simpson_rule(const integration_fn fn, const double start, const double end, const size_t segments) {
    double step = (end - start) / segments;
    double sum = fn(start) + fn(end);
    for (size_t i = 1; i < segments; i++) {
        sum += ((i % 2) ? 4.0 : 2.0) * fn(start + i * step);
    }
    return (step * sum) / 3.0;
}

double simpson_rule_38(const integration_fn fn, const double start, const double end, size_t segments) {
    // segments must be multiple of 3!
    segments = 3 * ((segments + 2) / 3);
    double step = (end - start) / segments;
    double sum = fn(start) + fn(end);
    for (size_t i = 1; i < segments; i++) {
        sum += ((i % 3) ? 3.0 : 2.0) * fn(start + i * step);
    }
    return (3.0 / 8.0) * (step * sum);
}

void gauss_legendre_coefficients(const size_t order, Vector *const t, Vector *const c) {
    if ((order == 0) || (t == NULL) || (c == NULL)) {
        return;  // Invalid operation
    }
    *t = polynomial_legendre_roots(order);
    if (!vector_is_valid(*t)) {
        return;
    }
    Matrix A = matrix_alloc(order, order);
    *c = vector_alloc(order);
    if (!matrix_is_valid(A) || !vector_is_valid(*c)) {
        matrix_dealloc(&A);
        vector_dealloc(c);
        return;
    }
    for (size_t i = 0; i < A.rows; i++) {
        for (size_t j = 0; j < A.cols; j++) {
            matrix_set(A, i, j, power(t->data[j], i));
        }
        c->data[i] = (i % 2) ? 0.0 : (2.0 / ((double)(i + 1)));
    }
    gaussian_elimination_over(A, *c);
    matrix_dealloc(&A);
}

double gauss_legendre_quadrature(const integration_fn fn, const double start, const double end, const size_t segments) {
    if (segments == 0) {
        return NAN;  // Invalid operation
    }
    Vector t = (Vector){0};
    Vector c = (Vector){0};
    gauss_legendre_coefficients(segments, &t, &c);
    double quadrature = NAN;
    if (vector_is_valid(t) && vector_is_valid(c)) {
        double sum = 0.0;
        for (size_t i = 0; i < t.len; i++) {
            const double xi = ((end - start) / 2.0) * t.data[i] + ((end + start) / 2.0);
            sum += c.data[i] * fn(xi);
        }
        quadrature = ((end - start) / 2.0) * sum;
    }
    vector_dealloc(&t);
    vector_dealloc(&c);
    return quadrature;
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