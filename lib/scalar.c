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

#include "scalar.h"

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>

#define MAX_ITERATIONS 10000
#define PRECISION 1e-10

void swap(double *const a, double *const b) {
    if ((a != NULL) && (b != NULL)) {
        const double temp = *a;
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

bool are_close(const double a, const double b, const double delta) {
    return (fabs(a - b) < delta);
}

double maximum(const double a, const double b) {
    return ((a > b) ? a : b);
}

double minimum(const double a, const double b) {
    return ((a < b) ? a : b);
}

double sign(const double value) {
    return ((value > 0) ? 1.0 : -1.0);
}

double random_number(const double min, const double max) {
    const double rand_unitary = ((double)rand()) / ((double)RAND_MAX);
    return rand_unitary * (max - min) + min;
}

// My own square root function, so I don't need to link with -lm,
// avoiding any dependencies
double square_root(const double value) {
    if (fabs(value) < PRECISION) {
        return 0.0;  // square_root(0) = 0
    } else if (value < 0) {
        return NAN;
    }
    double x = value;
    // Use the Newton-Raphson method to find the root
    // of the function f(x) = x^2 - value
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double delta = (value / x - x) / 2.0;
        x += delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return x;
}

// My own root function, so I don't need to link with -lm,
// avoiding any dependencies
double root(const double value, const uint64_t n) {
    if (n == 1) {
        return value;
    } else if ((n == 0) && (fabs(value) < PRECISION)) {
        return NAN;
    } else if (n == 0) {
        return 1.0;
    } else if (fabs(value) < PRECISION) {
        return 0.0;  // root(0) = 0
    } else if ((n % 2 == 0) && (value < 0)) {
        return NAN;
    }
    double x = value;
    // Use the Newton-Raphson method to find the root
    // of the function f(x) = x^n - value
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double delta = (value / power(x, (n - 1)) - x) / ((double)n);
        x += delta;
        if (fabs(delta) < PRECISION) {
            break;
        }
    }
    return x;
}

// My own power function, so I don't need to link with -lm,
// avoiding any dependencies
double power(const double base, uint64_t expoent) {
    double result = 1.0;
    while (expoent != 0) {
        result *= base;
        expoent--;
    }
    return result;
}

// My own exponential function, so I don't need to link with -lm,
// avoiding any dependencies
double exponential(const double value) {
    double term = 1.0;
    double result = term;
    for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
        term *= value / ((double)i);
        result += term;
        if (fabs(term) < PRECISION) {
            break;
        }
    }
    return result;
}

// My own sine function, so I don't need to link with -lm,
// avoiding any dependencies
double sine(double value) {
    value = value - 2 * PI * floor(value / (2 * PI));
    double term = value;
    double result = term;
    for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
        term *= -value * value / ((double)((2 * i + 1) * (2 * i)));
        result += term;
        if (fabs(term) < PRECISION) {
            break;
        }
    }
    return result;
}

// My own cosine function, so I don't need to link with -lm,
// avoiding any dependencies
double cosine(double value) {
    value = value - 2 * PI * floor(value / (2 * PI));
    double term = 1.0;
    double result = term;
    for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
        term *= -value * value / ((double)((2 * i) * (2 * i - 1)));
        result += term;
        if (fabs(term) < PRECISION) {
            break;
        }
    }
    return result;
}

// My own tangent function, so I don't need to link with -lm,
// avoiding any dependencies
double tangent(const double value) {
    return (sine(value) / cosine(value));
}

// My own arctangent function, so I don't need to link with -lm,
// avoiding any dependencies
double arctangent(const double value) {
    if (fabs(value) == 1.0) {
        return (sign(value) * PI / 4.0);
    } else if (fabs(value) < 1.0) {
        double result = value;
        for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
            const double term = (i % 2 ? -1.0 : 1.0) * power(value, 2 * i + 1) / ((double)(2 * i + 1));
            result += term;
            if (fabs(term) < PRECISION) {
                break;
            }
        }
        return result;
    } else {  // if fabs(value) >= 1.0
        double result = PI / 2.0 * (value > 1.0 ? 1.0 : -1.0);
        for (uint64_t i = 1; i < MAX_ITERATIONS; i++) {
            const double term = (i % 2 ? -1.0 : 1.0) / (power(value, 2 * i - 1) * ((double)(2 * i - 1)));
            result += term;
            if (fabs(term) < PRECISION) {
                break;
            }
        }
        return result;
    }
}

double atan2(const double y, const double x) {
    if ((x == 0.0) && (y == 0.0)) {
        return 0.0;
    } else if (x == 0.0) {
        return (sign(y) * PI / 2.0);
    } else if (x > 0.0) {
        return arctangent(y / x);
    } else if (y >= 0.0) {
        return (arctangent(y / x) + PI);
    } else {
        return (arctangent(y / x) - PI);
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