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

#include "complex-vector.h"
#include "complex.h"
#include "scalar.h"
#include "sorting.h"

#define MAX_ITERATIONS 10000
#define PRECISION 1e-12

void polynomial_print(const Vector polynomial, const char *const name, char x) {
    if (!isalpha(x)) {
        x = 'x';
    }
    if ((name != NULL) && (*name != '\0')) {
        printf("%s(%c) = ", name, x);
    }
    size_t printed_terms = 0;
    for (size_t i = (polynomial.len - 1); i < polynomial.len; i--) {
        if (are_close(polynomial.data[i], 0.0, PRECISION)) {
            continue;
        } else if (polynomial.data[i] < 0.0) {
            printf(" - ");
        } else if ((polynomial.data[i] > 0.0) && (printed_terms > 0)) {
            printf(" + ");
        }
        if (i == 0) {
            printf("%lg", fabs(polynomial.data[i]));
        } else {
            if (!are_close(fabs(polynomial.data[i]), 1.0, PRECISION)) {
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
Vector polynomial_mul(const Vector poly1, const Vector poly2) {
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
// Remember to free the 'result' vector after calling this function!
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
    if (!vector_is_valid(*result)) {
        return NAN;
    }
    result->data[result->len - 1] = polynomial.data[result->len];
    for (size_t i = (result->len - 2); i < result->len; i--) {
        result->data[i] = polynomial.data[i + 1] + r * result->data[i + 1];
    }
    // return the remainder
    return (polynomial.data[0] + r * result->data[0]);
}

// Computes the division of the polynomial by (x^2 + a*x + b)
// Remember to free the returned vector and 'result' after calling this function!
Vector polynomial_quadratic_division(const Vector polynomial, const double a, const double b, Vector *const result) {
    if (result == NULL) {
        return (Vector){0};  // Invalid operation
    }
    if (polynomial.len <= 2) {
        *result = (Vector){0};
        return vector_copy(polynomial);
    }
    *result = vector_alloc(polynomial.len - 2);
    if (!vector_is_valid(*result)) {
        return (Vector){0};
    }
    result->data[result->len - 1] = polynomial.data[result->len + 1];
    if (result->len > 2) {
        result->data[result->len - 2] = polynomial.data[result->len] - a * result->data[result->len - 1];
        for (size_t i = (result->len - 3); i < result->len; i--) {
            result->data[i] = polynomial.data[i + 2] - a * result->data[i + 1] - b * result->data[i + 2];
        }
    }
    // return the remainder
    Vector remainder = vector_alloc(2);
    if (vector_is_valid(remainder)) {
        remainder.data[1] = polynomial.data[1] - a * result->data[0] - ((result->len >= 2) ? (b * result->data[1]) : 0.0);
        remainder.data[0] = polynomial.data[0] - b * result->data[0];
    }
    return remainder;
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

// Polynomial differentiation using Ruffini's rule
double polynomial_diff(const Vector polynomial, const size_t order, const double x) {
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

// Remember to free the returned vector after calling this function!
Vector polynomial_ruffini_residuals(const Vector p, const double x) {
    Vector residuals = vector_init(p.len, 0.0);
    if (!vector_is_valid(residuals)) {
        return (Vector){0};
    }
    for (size_t i = p.len - 1; i < p.len; i--) {
        residuals.data[0] = p.data[i] + x * residuals.data[0];
        for (size_t j = 1; j <= i; j++) {
            residuals.data[j] = residuals.data[j - 1] + x * residuals.data[j];
        }
    }
    return residuals;
}

// Horner's method
Complex polynomial_complex_evaluation(const Vector polynomial, const Complex c) {
    Complex y = complex_init(0.0, 0.0);
    if (polynomial.len > 0) {
        for (size_t i = polynomial.len - 1; i < polynomial.len; i--) {
            y = complex_sum(complex_init(polynomial.data[i], 0.0), complex_mul(c, y));
        }
    }
    return y;
}

// Polynomial differentiation using Ruffini's rule
Complex polynomial_complex_first_diff(const Vector polynomial, const Complex c) {
    if (polynomial.len <= 1) {
        return complex_init(0.0, 0.0);
    }
    Complex first_div_remainder = complex_init(0.0, 0.0);
    Complex second_div_remainder = complex_init(0.0, 0.0);
    for (size_t i = polynomial.len - 1; i > 0; i--) {
        first_div_remainder = complex_sum(complex_init(polynomial.data[i], 0.0), complex_mul(c, first_div_remainder));
        second_div_remainder = complex_sum(first_div_remainder, complex_mul(c, second_div_remainder));
    }
    return second_div_remainder;
}

// Polynomial differentiation using Ruffini's rule
Complex polynomial_complex_diff(const Vector polynomial, const size_t order, const Complex c) {
    if (polynomial.len <= order) {
        return complex_init(0.0, 0.0);
    }
    Complex_Vector remainders = complex_vector_init(order + 1, complex_init(0.0, 0.0));
    if (!complex_vector_is_valid(remainders)) {
        return complex_init(NAN, NAN);
    }
    for (size_t i = polynomial.len - 1; (i + 1) > order; i--) {
        remainders.data[0] = complex_sum(complex_init(polynomial.data[i], 0.0), complex_mul(c, remainders.data[0]));
        for (size_t j = 1; j < remainders.len; j++) {
            remainders.data[j] = complex_sum(remainders.data[j - 1], complex_mul(c, remainders.data[j]));
        }
    }
    const Complex diff = complex_mul_scalar(factorial(order), remainders.data[remainders.len - 1]);
    complex_vector_dealloc(&remainders);
    return diff;
}

// Remember to free the returned vector after calling this function!
Complex_Vector polynomial_complex_ruffini_residuals(const Vector p, const Complex x) {
    Complex_Vector residuals = complex_vector_init(p.len, complex_init(0.0, 0.0));
    if (!complex_vector_is_valid(residuals)) {
        return (Complex_Vector){0};
    }
    for (size_t i = p.len - 1; i < p.len; i--) {
        residuals.data[0] = complex_sum(complex_init(p.data[i], 0.0), complex_mul(x, residuals.data[0]));
        for (size_t j = 1; ((j <= i) && (j < residuals.len)); j++) {
            residuals.data[j] = complex_sum(residuals.data[j - 1], complex_mul(x, residuals.data[j]));
        }
    }
    return residuals;
}

// Remember to free the returned vector after calling this function!
Vector polynomial_legendre(const size_t order) {
    if (order == 0) {
        Vector p0 = vector_alloc(1);
        if (vector_is_valid(p0)) {
            p0.data[0] = 1.0;
        }
        return p0;
    } else if (order == 1) {
        Vector p1 = vector_alloc(2);
        if (vector_is_valid(p1)) {
            p1.data[0] = 0.0;
            p1.data[1] = 1.0;
        }
        return p1;
    }
    const size_t num_of_p = 3;
    Vector p[num_of_p];
    for (size_t k = 0; k < num_of_p; k++) {
        p[k] = vector_init(order + 1, 0.0);
    }
    p[0].data[0] = 1.0;
    p[1].data[0] = 0.0;
    p[1].data[1] = 1.0;
    // P(n) = ((2n-1)/n)*x*P(n-1) - ((n-1)/n)*P(n-2)
    for (size_t k = 2; k <= order; k++) {
        const double a = (((double)(2 * k - 1) / ((double)k)));
        const double b = (((double)(k - 1)) / ((double)k));
        p[k % num_of_p].data[0] = -b * p[(k + 1) % num_of_p].data[0];
        for (size_t i = 1; i <= k; i++) {
            p[k % num_of_p].data[i] = a * p[(k + 2) % num_of_p].data[i - 1] - b * p[(k + 1) % num_of_p].data[i];
        }
    }
    vector_dealloc(&p[(order + 1) % num_of_p]);
    vector_dealloc(&p[(order + 2) % num_of_p]);
    return p[order % num_of_p];
}

Vector polynomial_legendre_roots(const size_t order) {
    Vector p = polynomial_legendre(order);
    if (!vector_is_valid(p)) {
        return (Vector){0};
    }
    Vector roots = polynomial_find_real_distinct_roots(p);
    vector_dealloc(&p);
    return roots;
}

// Remember to free the returned vector after calling this function!
Vector polynomial_chebyshev_first_kind(const size_t order) {
    if (order == 0) {
        Vector p0 = vector_alloc(1);
        if (vector_is_valid(p0)) {
            p0.data[0] = 1.0;
        }
        return p0;
    } else if (order == 1) {
        Vector p1 = vector_alloc(2);
        if (vector_is_valid(p1)) {
            p1.data[0] = 0.0;
            p1.data[1] = 1.0;
        }
        return p1;
    }
    const size_t num_of_p = 3;
    Vector p[num_of_p];
    for (size_t k = 0; k < num_of_p; k++) {
        p[k] = vector_init(order + 1, 0.0);
    }
    p[0].data[0] = 1.0;
    p[1].data[0] = 0.0;
    p[1].data[1] = 1.0;
    // T(n) = 2*x*T(n-1) - T(n-2)
    for (size_t k = 2; k <= order; k++) {
        p[k % num_of_p].data[0] = -p[(k + 1) % num_of_p].data[0];
        for (size_t i = 1; i <= k; i++) {
            p[k % num_of_p].data[i] = 2.0 * p[(k + 2) % num_of_p].data[i - 1] - p[(k + 1) % num_of_p].data[i];
        }
    }
    vector_dealloc(&p[(order + 1) % num_of_p]);
    vector_dealloc(&p[(order + 2) % num_of_p]);
    return p[order % num_of_p];
}

// Remember to free the returned vector after calling this function!
Vector polynomial_chebyshev_second_kind(const size_t order) {
    if (order == 0) {
        Vector p0 = vector_alloc(1);
        if (vector_is_valid(p0)) {
            p0.data[0] = 1.0;
        }
        return p0;
    } else if (order == 1) {
        Vector p1 = vector_alloc(2);
        if (vector_is_valid(p1)) {
            p1.data[0] = 0.0;
            p1.data[1] = 2.0;
        }
        return p1;
    }
    const size_t num_of_p = 3;
    Vector p[num_of_p];
    for (size_t k = 0; k < num_of_p; k++) {
        p[k] = vector_init(order + 1, 0.0);
    }
    p[0].data[0] = 1.0;
    p[1].data[0] = 0.0;
    p[1].data[1] = 2.0;
    // U(n) = 2*x*U(n-1) - U(n-2)
    for (size_t k = 2; k <= order; k++) {
        p[k % num_of_p].data[0] = -p[(k + 1) % num_of_p].data[0];
        for (size_t i = 1; i <= k; i++) {
            p[k % num_of_p].data[i] = 2.0 * p[(k + 2) % num_of_p].data[i - 1] - p[(k + 1) % num_of_p].data[i];
        }
    }
    vector_dealloc(&p[(order + 1) % num_of_p]);
    vector_dealloc(&p[(order + 2) % num_of_p]);
    return p[order % num_of_p];
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
    if (polynomial.len <= 1) {
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
}

// Implemeted according to the book Cáculo Numérico Computacional, pg. 267
// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_limit(const Vector p, size_t length) {
    // Computes the first limit
    double least_dif = INFINITY;
    double greater_coef = 0.0;
    for (size_t i = (p.len - length); i < p.len; i++) {
        greater_coef = maximum(greater_coef, fabs(p.data[i]));
        for (size_t j = (p.len - length); (j < p.len); j++) {
            if (i != j) {
                const double dif = fabs(fabs(p.data[i]) - fabs(p.data[j]));
                least_dif = minimum(least_dif, dif);
            }
        }
    }
    double limit1 = 0.1 * (least_dif / greater_coef);
    if (limit1 < PRECISION) {
        limit1 = 0.1;
    }
    // Computes the second limit
    int64_t decimals = 0;
    for (size_t i = (p.len - length); i < p.len; i++) {
        double coef = p.data[i];
        double dif = 1.0;
        int64_t cont = -1;
        while (dif > 1e-8) {
            dif = fabs(coef - round(coef));
            coef = 10.0 * coef;
            cont++;
        }
        decimals = (decimals > cont) ? decimals : cont;
    }
    const double limit2 = power(0.1, decimals);
    // The overall limit is
    const double limit = maximum(limit1 * limit2, 1e-7);
    return limit;
}

// This function only uses the first length coefficients of the polynomial in calculations
void polynomial_residuals_over(const Vector p, const Complex x, size_t length, Complex_Vector residuals) {
    complex_vector_init_over(residuals, complex_init(0.0, 0.0));
    for (size_t i = p.len - 1; ((i < p.len) && (i >= (p.len - length))); i--) {
        residuals.data[0] = complex_sum(complex_init(p.data[i], 0.0), complex_mul(x, residuals.data[0]));
        for (size_t j = 1; j <= (i - residuals.len + length); j++) {
            residuals.data[j] = complex_sum(residuals.data[j - 1], complex_mul(x, residuals.data[j]));
        }
    }
}

// Find root multiplicity
size_t polynomial_root_multiplicity(const Vector p, const Complex root) {
    Complex_Vector residuals = complex_vector_alloc(p.len);
    if (!complex_vector_is_valid(residuals)) {
        return 0;
    }
    polynomial_residuals_over(p, root, p.len, residuals);
    const double limit = polynomial_limit(p, p.len);
    size_t multiplicity = 1;
    double remainders = complex_modulus(residuals.data[0]) + complex_modulus(residuals.data[1]);
    while (remainders < limit) {
        multiplicity++;
        remainders += complex_modulus(residuals.data[multiplicity]);
    }
    complex_vector_dealloc(&residuals);
    return multiplicity;
}

// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_eval_over(const Vector polynomial, const double x, const size_t length) {
    double y = 0.0;
    for (size_t j = polynomial.len - 1; ((j >= polynomial.len - length) && (j < polynomial.len)); j--) {
        y = polynomial.data[j] + x * y;
    }
    return y;
}

// Computes the division of the polynomial p by (x-root) using Ruffini's rule
// This function only uses the first length coefficients of the polynomial in calculations
void polynomial_div_over(const Vector polynomial, const double r, const size_t length) {
    for (size_t j = (polynomial.len - 2); ((j >= polynomial.len - length) && (j < polynomial.len)); j--) {
        polynomial.data[j] = polynomial.data[j] + r * polynomial.data[j + 1];
    }
}

// Computes the division of the polynomial p by (x^2 + a*x + b)
// This function only uses the first length coefficients of the polynomial in calculations
void polynomial_quadratic_div_over(const Vector polynomial, const double a, const double b, const size_t length) {
    if (length > 2) {
        polynomial.data[polynomial.len - 2] = polynomial.data[polynomial.len - 2] - a * polynomial.data[polynomial.len - 1];
        for (size_t j = (polynomial.len - 3); ((j >= polynomial.len - length) && (j < polynomial.len)); j--) {
            polynomial.data[j] = polynomial.data[j] - a * polynomial.data[j + 1] - b * polynomial.data[j + 2];
        }
    }
}

// Cauchy's upper bound
// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_cauchy_upper_bound_over(const Vector polynomial, const size_t length) {
    double result = 0.0;
    for (size_t i = (polynomial.len - length); (i + 1) < polynomial.len; i++) {
        result = maximum(result, fabs(polynomial.data[i]));
    }
    return 1.0 + (result / fabs(polynomial.data[polynomial.len - 1]));
}

// Cauchy's lower bound
// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_cauchy_lower_bound_over(const Vector polynomial, const size_t length) {
    double result = 0.0;
    for (size_t i = (polynomial.len - length + 1); i < polynomial.len; i++) {
        result = maximum(result, fabs(polynomial.data[i]));
    }
    result = 1.0 + (result / fabs(polynomial.data[polynomial.len - length]));
    return (1.0 / result);
}

// Lagrange's upper bound
// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_lagrange_upper_bound_over(const Vector polynomial, const size_t length) {
    double result = 0.0;
    for (size_t i = (polynomial.len - length); (i + 1) < polynomial.len; i++) {
        result += fabs(polynomial.data[i]);
    }
    result /= fabs(polynomial.data[polynomial.len - 1]);
    return maximum(1.0, result);
}

// Lagrange's lower bound
// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_lagrange_lower_bound_over(const Vector polynomial, const size_t length) {
    double result = 0.0;
    for (size_t i = (polynomial.len - length + 1); i < polynomial.len; i++) {
        result += fabs(polynomial.data[i]);
    }
    result /= fabs(polynomial.data[polynomial.len - length]);
    return (1.0 / maximum(1.0, result));
}

// Cauchy's upper quota
// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_cauchy_upper_quota_over(const Vector polynomial, const size_t length) {
    double x = 0.0;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double previous_x = x;
        {
            x = 0.0;
            for (size_t i = 0; (i + 1) < length; i++) {
                x += fabs(polynomial.data[polynomial.len - length + i]) * power(previous_x, i);
            }
            x /= fabs(polynomial.data[polynomial.len - 1]);
            x = root(x, (length - 1));
        }
        if (fabs(x - previous_x) < PRECISION) {
            break;
        }
    }
    return x;
}

// Cauchy's lower quota
// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_cauchy_lower_quota_over(const Vector polynomial, const size_t length) {
    double x = 0.0;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        const double previous_x = x;
        {
            x = 0.0;
            for (size_t i = 1; i < length; i++) {
                x += fabs(polynomial.data[polynomial.len - length + i]) * power(previous_x, (length - i - 1));
            }
            x /= fabs(polynomial.data[polynomial.len - length]);
            x = root(x, (length - 1));
        }
        if (fabs(x - previous_x) < PRECISION) {
            break;
        }
    }
    return (1.0 / x);
}

// Kojima's upper bound
// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_kojima_upper_bound_over(const Vector polynomial, const size_t length) {
    double max = 0.0, second_max = 0.0;
    for (size_t i = (polynomial.len - length); (i + 1) < polynomial.len; i++) {
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
// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_kojima_lower_bound_over(const Vector polynomial, const size_t length) {
    double max = 0.0, second_max = 0.0;
    for (size_t i = (polynomial.len - length + 1); i < polynomial.len; i++) {
        const double value = root(fabs(polynomial.data[i] / polynomial.data[polynomial.len - length]), (i - polynomial.len + length));
        if (max < value) {
            second_max = max;
            max = value;
        } else if (second_max < value) {
            second_max = value;
        }
    }
    return (1.0 / (max + second_max));
}

// This function only uses the first length coefficients of the polynomial in calculations
void polynomial_root_bounds_over(const Vector polynomial, double *const min, double *const max, const size_t length) {
    if (max != NULL) {
        *max = polynomial_cauchy_upper_bound_over(polynomial, length);
        *max = minimum(*max, polynomial_lagrange_upper_bound_over(polynomial, length));
        *max = minimum(*max, polynomial_cauchy_upper_quota_over(polynomial, length));
        *max = minimum(*max, polynomial_kojima_upper_bound_over(polynomial, length));
    }
    if (min != NULL) {
        *min = polynomial_cauchy_lower_bound_over(polynomial, length);
        *min = maximum(*min, polynomial_lagrange_lower_bound_over(polynomial, length));
        *min = maximum(*min, polynomial_cauchy_lower_quota_over(polynomial, length));
        *min = maximum(*min, polynomial_kojima_lower_bound_over(polynomial, length));
    }
}

// This function only uses the first length coefficients of the polynomial in calculations
Complex polynomial_root_guess_over(const Vector polynomial, const size_t length) {
    double min, max;
    polynomial_root_bounds_over(polynomial, &min, &max, length);
    if (are_close((min - max), 0.0, PRECISION)) {
        min *= 0.99;
        max *= 1.01;
    }
    //  Try to find a real guess for the root
    const double step = (max - min) / ((double)(10 * length));
    double a = -max - PRECISION;
    double b = a + step;
    while (b < (PRECISION - min)) {
        const double fa = polynomial_eval_over(polynomial, a, length);
        const double fb = polynomial_eval_over(polynomial, b, length);
        if (fa * fb <= 0.0) {
            return complex_init((a + b) / 2.0, 0.0);
        }
        a = b;
        b = b + step;
    }
    a = min - PRECISION;
    b = a + step;
    while (b < (PRECISION + max)) {
        const double fa = polynomial_eval_over(polynomial, a, length);
        const double fb = polynomial_eval_over(polynomial, b, length);
        if (fa * fb <= 0.0) {
            return complex_init((a + b) / 2.0, 0.0);
        }
        a = b;
        b = b + step;
    }
    // If the polynomial has odd order, there is at least one real root
    if ((length % 2) == 0) {
        return complex_init(square_root(max * min), 0.0);
    }
    // If didn't found a real root guess, try a complex initial guess
    const double guess = square_root((max * min) / 2.0);
    return complex_init(guess, guess);
}

// This function only uses the first length coefficients of the polynomial in calculations
size_t polynomial_find_root_over(const Vector p, Complex *const root, size_t length, Complex_Vector residuals) {
    if (root == NULL) {
        return 0;  // Invalid operation
    }
    const double limit = polynomial_limit(p, length);
    size_t multiplicity;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        // Find root multiplicity
        multiplicity = 1;
        polynomial_residuals_over(p, *root, length, residuals);
        double remainders = complex_modulus(residuals.data[0]) + complex_modulus(residuals.data[1]);
        while (remainders < limit) {
            multiplicity++;
            remainders += complex_modulus(residuals.data[multiplicity]);
        }
        // Use Newton method for root finding - Birge-Vieta method
        const Complex delta = complex_div(residuals.data[multiplicity - 1], complex_mul_scalar(multiplicity, residuals.data[multiplicity]));
        *root = complex_sub(*root, delta);
        if ((complex_modulus(residuals.data[0]) + complex_modulus(delta)) < PRECISION) {
            break;
        }
    }
    return multiplicity;
}

void polynomial_root_refinement_over(const Vector p, Complex *const root, const size_t multiplicity, Complex_Vector residuals) {
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        polynomial_residuals_over(p, *root, p.len, residuals);
        // Use Newton method for root finding - Birge-Vieta method
        const Complex delta = complex_div(residuals.data[multiplicity - 1], complex_mul_scalar(multiplicity, residuals.data[multiplicity]));
        *root = complex_sub(*root, delta);
        if ((complex_modulus(residuals.data[0]) + complex_modulus(delta)) < PRECISION) {
            break;
        }
    }
}

// Remember to free the returned vector after calling this function!
Complex_Vector polynomial_find_roots(const Vector polynomial) {
    if (polynomial.len <= 1) {
        return (Complex_Vector){0};  // Invalid operation
    }
    Vector p = vector_alloc(polynomial.len);
    Complex_Vector roots = complex_vector_alloc(polynomial.len - 1);
    Complex_Vector residuals = complex_vector_alloc(polynomial.len);
    if ((!vector_is_valid(p)) || (!complex_vector_is_valid(roots)) || (!complex_vector_is_valid(residuals))) {
        vector_dealloc(&p);
        complex_vector_dealloc(&roots);
        complex_vector_dealloc(&residuals);
        return (Complex_Vector){0};
    }
    // Normalization of the coefficients of the polynomial
    p.data[polynomial.len - 1] = 1.0;
    for (size_t i = 0; (i + 1) < polynomial.len; i++) {
        p.data[i] = polynomial.data[i] / polynomial.data[polynomial.len - 1];
    }
    /*
    OBS: In order to avoid realocating memory over and over again, we store
         the result of the polynomial division of p by (x - root) in p as well.
         The length of the p polynomial is (p.len - i).
    */
    // Check for null roots
    size_t root_index = 0;
    while (are_close(p.data[root_index], 0.0, PRECISION)) {
        roots.data[root_index] = complex_init(p.data[root_index], 0.0);
        root_index++;
    }
    // Find remaining roots
    while (root_index < roots.len) {
        const size_t length = (p.len - root_index);
        Complex root = polynomial_root_guess_over(p, length);
        size_t multiplicity = polynomial_find_root_over(p, &root, length, residuals);
        polynomial_root_refinement_over(polynomial, &root, multiplicity, residuals);
        while ((multiplicity--) && (root_index < roots.len)) {
            if (complex_is_real(root)) {
                //  Computes the division of the polynomial p by (x-root) using Ruffini's rule
                polynomial_div_over(p, root.real, length);
                // Add the root found to the result vector
                roots.data[root_index] = root;
                root_index++;
            } else {
                // Computes the division of the polynomial p by (x^2 + a*x + b)
                const double a = -2.0 * root.real;
                const double b = root.real * root.real + root.imag * root.imag;
                polynomial_quadratic_div_over(p, a, b, length);
                // Add the roots found to the result vector
                roots.data[root_index] = root;
                root_index++;
                if (root_index >= roots.len) {
                    break;
                }
                roots.data[root_index] = complex_conjugate(root);
                root_index++;
            }
        }
    }
    complex_vector_dealloc(&residuals);
    vector_dealloc(&p);
    return roots;
}

// This function only uses the first length coefficients of the polynomial in calculations
void polynomial_real_residuals_over(const Vector p, const double x, size_t length, Vector residuals) {
    vector_init_over(residuals, 0.0);
    for (size_t i = p.len - 1; ((i < p.len) && (i >= (p.len - length))); i--) {
        residuals.data[0] = p.data[i] + x * residuals.data[0];
        for (size_t j = 1; j <= (i - residuals.len + length); j++) {
            residuals.data[j] = residuals.data[j - 1] + x * residuals.data[j];
        }
    }
}

// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_real_root_guess_over(const Vector polynomial, const size_t length) {
    double min, max;
    polynomial_root_bounds_over(polynomial, &min, &max, length);
    if (are_close((min - max), 0.0, PRECISION)) {
        min *= 0.99;
        max *= 1.01;
    }
    //  Try to find a real guess for the root
    const double step = (max - min) / ((double)(10 * length));
    double a = -max - PRECISION;
    double b = a + step;
    while (b < (PRECISION - min)) {
        const double fa = polynomial_eval_over(polynomial, a, length);
        const double fb = polynomial_eval_over(polynomial, b, length);
        if (fa * fb <= 0.0) {
            return ((a + b) / 2.0);
        }
        a = b;
        b = b + step;
    }
    a = min - PRECISION;
    b = a + step;
    while (b < (PRECISION + max)) {
        const double fa = polynomial_eval_over(polynomial, a, length);
        const double fb = polynomial_eval_over(polynomial, b, length);
        if (fa * fb <= 0.0) {
            return ((a + b) / 2.0);
        }
        a = b;
        b = b + step;
    }
    return square_root(max * min);
}

// This function only uses the first length coefficients of the polynomial in calculations
size_t polynomial_find_real_root_over(const Vector p, double *const root, size_t length, Vector residuals) {
    if (root == NULL) {
        return 0;  // Invalid operation
    }
    const double limit = polynomial_limit(p, length);
    size_t multiplicity;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        // Find root multiplicity
        multiplicity = 1;
        polynomial_real_residuals_over(p, *root, length, residuals);
        double remainders = fabs(residuals.data[0]) + fabs(residuals.data[1]);
        while (remainders < limit) {
            multiplicity++;
            remainders += fabs(residuals.data[multiplicity]);
        }
        // Use Newton method for root finding - Birge-Vieta method
        const double delta = residuals.data[multiplicity - 1] / (((double)multiplicity) * residuals.data[multiplicity]);
        *root -= delta;
        if ((fabs(residuals.data[0]) + fabs(delta)) < PRECISION) {
            break;
        }
    }
    return multiplicity;
}

void polynomial_real_root_refinement_over(const Vector p, double *const root, const size_t multiplicity, Vector residuals) {
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        polynomial_real_residuals_over(p, *root, p.len, residuals);
        // Use Newton method for root finding - Birge-Vieta method
        const double delta = residuals.data[multiplicity - 1] / (((double)multiplicity) * residuals.data[multiplicity]);
        *root -= delta;
        if ((fabs(residuals.data[0]) + fabs(delta)) < PRECISION) {
            break;
        }
    }
}

// Remember to free the returned vector after calling this function!
// This function only works for polynomials which all roots are real
Vector polynomial_find_real_roots(const Vector polynomial) {
    if (polynomial.len <= 1) {
        return (Vector){0};  // Invalid operation
    }
    Vector p = vector_alloc(polynomial.len);
    Vector roots = vector_alloc(polynomial.len - 1);
    Vector residuals = vector_alloc(polynomial.len);
    if ((!vector_is_valid(p)) || (!vector_is_valid(roots)) || (!vector_is_valid(residuals))) {
        vector_dealloc(&p);
        vector_dealloc(&roots);
        vector_dealloc(&residuals);
        return (Vector){0};
    }
    // Normalization of the coefficients of the polynomial
    p.data[polynomial.len - 1] = 1.0;
    for (size_t i = 0; (i + 1) < polynomial.len; i++) {
        p.data[i] = polynomial.data[i] / polynomial.data[polynomial.len - 1];
    }
    /*
    OBS: In order to avoid realocating memory over and over again, we store
         the result of the polynomial division of p by (x - root) in p as well.
         The length of the p polynomial is (p.len - i).
    */
    // Check for null roots
    size_t root_index = 0;
    while (are_close(p.data[root_index], 0.0, PRECISION)) {
        roots.data[root_index] = p.data[root_index];
        root_index++;
    }
    // Find remaining roots
    while (root_index < roots.len) {
        const size_t length = (p.len - root_index);
        double root = polynomial_real_root_guess_over(p, length);
        size_t multiplicity = polynomial_find_real_root_over(p, &root, length, residuals);
        polynomial_real_root_refinement_over(polynomial, &root, multiplicity, residuals);
        while ((multiplicity--) && (root_index < roots.len)) {
            //  Computes the division of the polynomial p by (x-root) using Ruffini's rule
            polynomial_div_over(p, root, length);
            // Add the root found to the result vector
            roots.data[root_index] = root;
            root_index++;
        }
    }
    vector_dealloc(&residuals);
    vector_dealloc(&p);
    quicksort(roots);
    return roots;
}

// This function only uses the first length coefficients of the polynomial in calculations
double polynomial_first_diff_over(const Vector p, const double x, double *const y, size_t length) {
    double first_div_remainder = 0.0;
    double second_div_remainder = 0.0;
    for (size_t i = p.len - 1; ((i < p.len) && (i > (p.len - length))); i--) {
        first_div_remainder = p.data[i] + x * first_div_remainder;
        second_div_remainder = first_div_remainder + x * second_div_remainder;
    }
    if (y != NULL) {
        first_div_remainder = p.data[p.len - length] + x * first_div_remainder;
        *y = first_div_remainder;
    }
    return second_div_remainder;
}

// This function only uses the first length coefficients of the polynomial in calculations
void polynomial_find_real_root_newton_over(const Vector p, double *const root, size_t length) {
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        double y;
        const double diff = polynomial_first_diff_over(p, *root, &y, length);
        // Use Newton method for root finding - Birge-Vieta method
        const double delta = y / diff;
        *root -= delta;
        if ((fabs(y) + fabs(delta)) < PRECISION) {
            break;
        }
    }
}

// Remember to free the returned vector after calling this function!
// This function only works for polynomials which all roots are real
// and non-repetitive (distinct)
Vector polynomial_find_real_distinct_roots(const Vector polynomial) {
    if (polynomial.len <= 1) {
        return (Vector){0};  // Invalid operation
    }
    Vector p = vector_alloc(polynomial.len);
    Vector roots = vector_alloc(polynomial.len - 1);
    if ((!vector_is_valid(p)) || (!vector_is_valid(roots))) {
        vector_dealloc(&p);
        vector_dealloc(&roots);
        return (Vector){0};
    }
    // Normalization of the coefficients of the polynomial
    p.data[polynomial.len - 1] = 1.0;
    for (size_t i = 0; (i + 1) < polynomial.len; i++) {
        p.data[i] = polynomial.data[i] / polynomial.data[polynomial.len - 1];
    }
    /*
    OBS: In order to avoid realocating memory over and over again, we store
         the result of the polynomial division of p by (x - root) in p as well.
         The length of the p polynomial is (p.len - i).
    */
    // Check for null roots
    size_t root_index = 0;
    while (are_close(p.data[root_index], 0.0, PRECISION)) {
        roots.data[root_index] = p.data[root_index];
        root_index++;
    }
    // Find remaining roots
    while (root_index < roots.len) {
        const size_t length = (p.len - root_index);
        double root = polynomial_real_root_guess_over(p, length);
        polynomial_find_real_root_newton_over(p, &root, length);
        // Root refinement
        polynomial_find_real_root_newton_over(polynomial, &root, polynomial.len);
        //  Computes the division of the polynomial p by (x-root) using Ruffini's rule
        polynomial_div_over(p, root, length);
        // Add the root found to the result vector
        roots.data[root_index] = root;
        root_index++;
    }
    vector_dealloc(&p);
    quicksort(roots);
    return roots;
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