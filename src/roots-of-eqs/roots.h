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

#ifndef __ROOTS_H
#define __ROOTS_H

#include <math.h>

#define MAX_ITERATIONS 10000
#define PRECISION 1e-10

typedef double (*root_fn)(const double);

double bisection_method(const root_fn fn, double start, double end);
double fakepos_method(const root_fn fn, double start, double end);
double newton_raphson_method(const root_fn fn, const root_fn dfn, const double initial);
double secant_method(const root_fn fn, const double start, const double end);

#endif  // __ROOTS_H

//------------------------------------------------------------------------------
// IMPLEMENTATION
//------------------------------------------------------------------------------

#ifdef ROOTS_IMPLEMENTATION

double bisection_method(const root_fn fn, double start, double end) {
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double middle = (start + end) / 2.0;
        if (fabs(end - start) < PRECISION) {
            return middle;
        }
        if (fn(start) * fn(middle) < 0) {
            end = middle;
        } else {
            start = middle;
        }
    }
    return 0.0;
}

double fakepos_method(const root_fn fn, double start, double end) {
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double x = end - ((fn(end) * (end - start)) / (fn(end) - fn(start)));
        if (fabs(end - start) < PRECISION) {
            return x;
        }
        if (fn(start) * fn(end) < 0) {
            end = x;
        } else {
            start = x;
        }
    }
    return 0.0;
}

double newton_raphson_method(const root_fn fn, const root_fn dfn, const double initial) {
    double x = initial;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double delta = fn(x) / dfn(x);
        x -= delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return x;
}

double secant_method(const root_fn fn, const double start, const double end) {
    double previous_x = start;
    double x = end;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double delta = (fn(x) * (x - previous_x)) / (fn(x) - fn(previous_x));
        previous_x = x;
        x -= delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return x;
}

#endif  // ROOTS_IMPLEMENTATION

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