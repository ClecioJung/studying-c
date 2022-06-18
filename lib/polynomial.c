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

#include "polynomial.h"

#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>

#include "scalar.h"
#include "sorting.h"

#define MAX_ITERATIONS 10000
#define PRECISION 1e-10

void polynomial_print(const Vector polynomial, const char *const name, char x) {
    if (!isalpha(x)) {
        x = 'x';
    }
    if ((name != NULL) && (*name != '\0')) {
        printf("%s(%c) = ", name, x);
    }
    size_t i = polynomial.len;
    size_t printed_terms = 0;
    while (i > 0) {
        i--;
        if (polynomial.data[i] == 0.0) {
            continue;
        } else if (polynomial.data[i] < 0.0) {
            printf(" - ");
        } else if ((polynomial.data[i] > 0.0) && (printed_terms > 0)) {
            printf(" + ");
        }
        if (i == 0) {
            printf("%lg", fabs(polynomial.data[i]));
        } else {
            if (polynomial.data[i] != 1.0) {
                printf("%lg*", fabs(polynomial.data[i]));
            }
            printf("%c", x);
            if (i != 1) {
                printf("^%ld", i);
            }
        }
        printed_terms++;
    }
    if (printed_terms == 0) {
        printf("0");
    }
    printf("\n");
}

// Remember to free the returned vector after calling this function!
Vector polynomial_sum(const Vector poly1, const Vector poly2) {
    Vector bigger_poly, small_poly;
    if (poly1.len > poly2.len) {
        bigger_poly = poly1;
        small_poly = poly2;
    } else {
        bigger_poly = poly2;
        small_poly = poly1;
    }
    Vector result = vector_copy(bigger_poly);
    if (vector_is_valid(result)) {
        for (size_t i = 0; i < small_poly.len; i++) {
            result.data[i] += small_poly.data[i];
        }
    }
    return result;
}

// Remember to free the returned vector after calling this function!
Vector polynomial_multiply(const Vector poly1, const Vector poly2) {
    Vector result = vector_init(poly1.len + poly2.len - 1, 0.0);
    if (vector_is_valid(result)) {
        for (size_t i = 0; i < poly1.len; i++) {
            for (size_t j = 0; j < poly2.len; j++) {
                result.data[i + j] += poly1.data[i] * poly2.data[j];
            }
        }
    }
    return result;
}

bool polynomial_are_equal(const Vector poly1, const Vector poly2) {
    return vector_are_equal(poly1, poly2);
}

double polynomial_evaluation(const Vector polynomial, const double x) {
    double y = 0.0;
    for (size_t i = 0; i < polynomial.len; i++) {
        y += polynomial.data[i] * power(x, i);
    }
    return y;
}

// Horner's method
double polynomial_horner_evaluation(const Vector polynomial, const double x) {
    if (polynomial.len == 0) {
        return 0.0;
    }
    double y = 0.0;
    for (size_t i = polynomial.len - 1; i < polynomial.len; i--) {
        y = polynomial.data[i] + x * y;
    }
    return y;
}

// Ruffini's rule
// Computes the division of the polynomial by a binomial of the form (x – r)
// Remember to free the returned vector after calling this function!
double polynomial_ruffini_division(const Vector polynomial, const double r, Vector *const result) {
    if (result == NULL) {
        return NAN;  // Invalid operation
    }
    if (polynomial.len == 0) {
        *result = (Vector){0};
        return 0.0;
    }
    if (polynomial.len == 1) {
        *result = (Vector){0};
        return polynomial.data[0];
    }
    *result = vector_alloc(polynomial.len - 1);
    if (result == NULL) {
        return NAN;
    }
    size_t i = result->len - 1;
    result->data[i] = polynomial.data[i + 1];
    while (i > 0) {
        i--;
        result->data[i] = polynomial.data[i + 1] + r * result->data[i + 1];
    }
    // return the remainder;
    return (polynomial.data[0] + r * result->data[0]);
}

// Polynomial differentiation using Ruffini's rule
double polynomial_first_diff(const Vector polynomial, const double x) {
    if (polynomial.len <= 1) {
        return 0.0;
    }
    double first_div_remainder = 0.0;
    double second_div_remainder = 0.0;
    for (size_t i = polynomial.len - 1; i > 0; i--) {
        first_div_remainder = polynomial.data[i] + x * first_div_remainder;
        second_div_remainder = first_div_remainder + x * second_div_remainder;
    }
    return second_div_remainder;
}

double polynomial_diff(const Vector polynomial, const uint16_t order, const double x) {
    if (polynomial.len <= order) {
        return 0.0;
    }
    Vector remainders = vector_init(order + 1, 0.0);
    if (!vector_is_valid(remainders)) {
        return NAN;
    }
    for (size_t i = polynomial.len - 1; (i + 1) > order; i--) {
        remainders.data[0] = polynomial.data[i] + x * remainders.data[0];
        for (size_t j = 1; j < remainders.len; j++) {
            remainders.data[j] = remainders.data[j - 1] + x * remainders.data[j];
        }
    }
    const double diff = factorial(order) * remainders.data[remainders.len - 1];
    vector_dealloc(&remainders);
    return diff;
}

// Cauchy's upper bound
double polynomial_cauchy_upper_bound(const Vector polynomial) {
    if (polynomial.len == 0) {
        return NAN;  // Invalid operation
    }
    double result = 0.0;
    for (size_t i = 0; (i + 1) < polynomial.len; i++) {
        result = maximum(result, fabs(polynomial.data[i]));
    }
    return 1.0 + (result / fabs(polynomial.data[polynomial.len - 1]));
}

// Cauchy's lower bound
double polynomial_cauchy_lower_bound(const Vector polynomial) {
    if (polynomial.len == 0) {
        return NAN;  // Invalid operation
    }
    double result = 0.0;
    for (size_t i = 1; i < polynomial.len; i++) {
        result = maximum(result, fabs(polynomial.data[i]));
    }
    result = 1.0 + (result / fabs(polynomial.data[0]));
    return (1.0 / result);
}

// Lagrange's upper bound
double polynomial_lagrange_upper_bound(const Vector polynomial) {
    if (polynomial.len == 0) {
        return NAN;  // Invalid operation
    }
    double result = 0.0;
    for (size_t i = 0; (i + 1) < polynomial.len; i++) {
        result += fabs(polynomial.data[i]);
    }
    result /= fabs(polynomial.data[polynomial.len - 1]);
    return maximum(1.0, result);
}

// Lagrange's lower bound
double polynomial_lagrange_lower_bound(const Vector polynomial) {
    if (polynomial.len == 0) {
        return NAN;  // Invalid operation
    }
    double result = 0.0;
    for (size_t i = 1; i < polynomial.len; i++) {
        result += fabs(polynomial.data[i]);
    }
    result /= fabs(polynomial.data[0]);
    return (1.0 / maximum(1.0, result));
}

// Cauchy's upper quota
double polynomial_cauchy_upper_quota(const Vector polynomial) {
    if (polynomial.len == 0) {
        return NAN;  // Invalid operation
    }
    double x = 0.0;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double previous_x = x;
        {
            x = 0.0;
            for (size_t i = 0; (i + 1) < polynomial.len; i++) {
                x += fabs(polynomial.data[i]) * power(previous_x, i);
            }
            x /= fabs(polynomial.data[polynomial.len - 1]);
            x = root(x, (polynomial.len - 1));
        }
        if (fabs(x - previous_x) < PRECISION) {
            break;
        }
    }
    return x;
}

// Cauchy's lower quota
double polynomial_cauchy_lower_quota(const Vector polynomial) {
    if (polynomial.len == 0) {
        return NAN;  // Invalid operation
    }
    double x = 0.0;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double previous_x = x;
        {
            x = 0.0;
            for (size_t i = 1; i < polynomial.len; i++) {
                x += fabs(polynomial.data[i]) * power(previous_x, (polynomial.len - i - 1));
            }
            x /= fabs(polynomial.data[0]);
            x = root(x, (polynomial.len - 1));
        }
        if (fabs(x - previous_x) < PRECISION) {
            break;
        }
    }
    return (1.0 / x);
}

// Kojima's upper bound
double polynomial_kojima_upper_bound(const Vector polynomial) {
    if (polynomial.len == 0) {
        return NAN;  // Invalid operation
    }
    double max = 0.0, second_max = 0.0;
    for (size_t i = 0; (i + 1) < polynomial.len; i++) {
        const double value = root(fabs(polynomial.data[i] / polynomial.data[polynomial.len - 1]), polynomial.len - i - 1);
        if (max < value) {
            second_max = max;
            max = value;
        } else if (second_max < value) {
            second_max = value;
        }
    }
    return (max + second_max);
}

// Kojima's lower bound
double polynomial_kojima_lower_bound(const Vector polynomial) {
    if (polynomial.len == 0) {
        return NAN;  // Invalid operation
    }
    double max = 0.0, second_max = 0.0;
    for (size_t i = 1; i < polynomial.len; i++) {
        const double value = root(fabs(polynomial.data[i] / polynomial.data[0]), i);
        if (max < value) {
            second_max = max;
            max = value;
        } else if (second_max < value) {
            second_max = value;
        }
    }
    return (1.0 / (max + second_max));
}

void polynomial_root_bounds(const Vector polynomial, double *const min, double *const max) {
    if (polynomial.len == 0) {
        return;  // Invalid operation
    }
    if (max != NULL) {
        *max = polynomial_cauchy_upper_bound(polynomial);
        *max = minimum(*max, polynomial_lagrange_upper_bound(polynomial));
        *max = minimum(*max, polynomial_cauchy_upper_quota(polynomial));
        *max = minimum(*max, polynomial_kojima_upper_bound(polynomial));
    }
    if (min != NULL) {
        *min = polynomial_cauchy_lower_bound(polynomial);
        *min = maximum(*min, polynomial_lagrange_lower_bound(polynomial));
        *min = maximum(*min, polynomial_cauchy_lower_quota(polynomial));
        *min = maximum(*min, polynomial_kojima_lower_bound(polynomial));
    }
    if (fabs(*min - *max) < PRECISION) {
        *min *= 0.99;
        *max *= 1.01;
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