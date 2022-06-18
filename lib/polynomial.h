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

#ifndef __POLYNOMIAL_H
#define __POLYNOMIAL_H

#include <stdbool.h>
#include <stdint.h>

#include "vector.h"

void polynomial_print(const Vector polynomial, const char *const name, char x);
Vector polynomial_sum(const Vector poly1, const Vector poly2);
Vector polynomial_multiply(const Vector poly1, const Vector poly2);
bool polynomial_are_equal(const Vector poly1, const Vector poly2);
double polynomial_evaluation(const Vector polynomial, const double x);
double polynomial_horner_evaluation(const Vector polynomial, const double x);
double polynomial_ruffini_division(const Vector polynomial, const double r, Vector *const result);
double polynomial_first_diff(const Vector polynomial, const double x);
double polynomial_diff(const Vector polynomial, uint16_t order, const double x);
double polynomial_cauchy_upper_bound(const Vector polynomial);
double polynomial_cauchy_lower_bound(const Vector polynomial);
double polynomial_lagrange_upper_bound(const Vector polynomial);
double polynomial_lagrange_lower_bound(const Vector polynomial);
double polynomial_cauchy_upper_quota(const Vector polynomial);
double polynomial_cauchy_lower_quota(const Vector polynomial);
double polynomial_kojima_upper_bound(const Vector polynomial);
double polynomial_kojima_lower_bound(const Vector polynomial);
void polynomial_root_bounds(const Vector polynomial, double *const min, double *const max);

#endif  // __POLYNOMIAL_H

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