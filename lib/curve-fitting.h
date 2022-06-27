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

#ifndef __CURVE_FITTING_H
#define __CURVE_FITTING_H

#include "vector.h"

#define SPLINE_ORDER 3

typedef struct {
    double *x;
    double (*coefficients)[SPLINE_ORDER + 1];
    size_t len;
} Spline;

typedef enum {
    Spline_Clamped = 0,
    Spline_Quadratic,
    Spline_Natural,
} Spline_Type;

double polynomial_interpolation(const Vector x, const Vector y, const double value);
double lagrange_interpolation(const Vector x, const Vector y, const double value);
double gregory_newton_interpolation(const Vector x, const Vector y, const double value);
double linear_regression(const Vector x, const Vector y, double *const a, double *const b);
Vector polynomial_regression(const Vector x, const Vector y, const size_t order);
double r_squared(const Vector x, const Vector y, const Vector polynomial);

Vector polynomial_interpolation_vector(const Vector x, const Vector y, const Vector values);
Vector lagrange_interpolation_vector(const Vector x, const Vector y, const Vector values);
Vector gregory_newton_interpolation_vector(const Vector x, const Vector y, const Vector values);

Spline spline_alloc(const size_t len);
void spline_dealloc(Spline *const spline);
Spline spline_interpolation(const Vector x, const Vector y, const Spline_Type spline_type, const double S_start, const double S_end);
void spline_print(const Spline spline);
double spline_evaluation(const Spline spline, const double value);

#endif  // __CURVE_FITTING_H

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