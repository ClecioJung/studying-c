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
    for (size_t i = 2; i < residuals.len; i++) {
        residuals.data[i] *= (double)factorial(i);
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
Complex polynomial_complex_diff(const Vector polynomial, const uint16_t order, const Complex c) {
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
    for (size_t i = 2; i < residuals.len; i++) {
        residuals.data[i] = complex_mul_scalar(factorial(i), residuals.data[i]);
    }
    return residuals;
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

// Implemeted according to the book Cáculo Numérico Computacional, pg. 267
double polynomial_limit1(const Vector p) {
    double least_dif = INFINITY;
    double greater_coef = 0.0;
    for (size_t i = 0; i < p.len; i++) {
        greater_coef = maximum(greater_coef, fabs(p.data[i]));
        for (size_t j = 0; j < p.len; j++) {
            if (i != j) {
                const double dif = fabs(fabs(p.data[i]) - fabs(p.data[j]));
                least_dif = minimum(least_dif, dif);
            }
        }
    }
    double limit = 0.1 * (least_dif / greater_coef);
    if (limit < PRECISION) {
        limit = 0.1;
    }
    return limit;
}

// Implemeted according to the book Cáculo Numérico Computacional, pg. 267
double polynomial_limit2(const Vector p) {
    uint16_t decimals = 0;
    for (size_t i = 0; i < p.len; i++) {
        double coef = p.data[i];
        double dif = 1.0;
        int16_t cont = -1;
        while (dif > 1e-8) {
            dif = fabs(coef - round(coef));
            coef = 10.0 * coef;
            cont++;
        }
        decimals = (decimals > cont) ? decimals : cont;
    }
    const double limit = power(0.1, decimals);
    return limit;
}

// Implemeted according to the book Cáculo Numérico Computacional, pg. 267
double polynomial_limit(const Vector p) {
    const double limit1 = polynomial_limit1(p);
    const double limit2 = polynomial_limit2(p);
    const double limit = maximum(limit1 * limit2, 1e-7);
    return limit;
}

uint16_t polynomial_root_multiplicity(const Vector p, const Complex root) {
    const double limit = polynomial_limit(p);
    uint16_t multiplicity = 1;
    double remainders = complex_modulus(polynomial_complex_diff(p, 0, root)) +
                        complex_modulus(polynomial_complex_diff(p, 1, root));
    while (remainders < limit) {
        multiplicity++;
        remainders += complex_modulus(polynomial_complex_diff(p, multiplicity, root));
    }
    return multiplicity;
}

Complex polynomial_root_guess(const Vector polynomial) {
    double min, max;
    polynomial_root_bounds(polynomial, &min, &max);
    // Try to find a real guess for the root
    const double step = (max - min) / ((double)(10 * polynomial.len));
    double a = -max - PRECISION;
    double b = a + step;
    while (b < (PRECISION - min)) {
        const double fa = polynomial_horner_evaluation(polynomial, a);
        const double fb = polynomial_horner_evaluation(polynomial, b);
        if (fa * fb <= 0.0) {
            return complex_init((a + b) / 2.0, 0.0);
        }
        a = b;
        b = b + step;
    }
    a = min - PRECISION;
    b = a + step;
    while (b < (PRECISION + max)) {
        const double fa = polynomial_horner_evaluation(polynomial, a);
        const double fb = polynomial_horner_evaluation(polynomial, b);
        if (fa * fb <= 0.0) {
            return complex_init((a + b) / 2.0, 0.0);
        }
        a = b;
        b = b + step;
    }
    // If didn't found a real root guess, try a complex initial guess
    const double guess = square_root((max * min) / 2.0);
    return complex_init(guess, guess);
}

uint16_t polynomial_find_root(const Vector p, Complex *const root) {
    if (root == NULL) {
        return 0;  // Invalid operation
    }
    const double limit = polynomial_limit(p);
    uint16_t multiplicity;
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        // Find root multiplicity
        multiplicity = 1;
        Complex numerator = polynomial_complex_evaluation(p, *root);
        Complex denominator = polynomial_complex_first_diff(p, *root);
        double remainders = complex_modulus(numerator) + complex_modulus(denominator);
        while (remainders < limit) {
            multiplicity++;
            numerator = denominator;
            denominator = polynomial_complex_diff(p, multiplicity, *root);
            remainders += complex_modulus(denominator);
        }
        // Use Newton method for root finding - Birge-Vieta method
        const Complex delta = complex_div(numerator, complex_mul_scalar(multiplicity, denominator));
        *root = complex_sub(*root, delta);
        if (complex_modulus(delta) < PRECISION) {
            break;
        }
    }
    return multiplicity;
}

void polynomial_root_refinement(const Vector p, Complex *const root, const uint16_t multiplicity) {
    for (size_t k = 0; k < MAX_ITERATIONS; k++) {
        // Use Newton method for root finding - Birge-Vieta method
        Complex numerator = polynomial_complex_diff(p, (multiplicity - 1), *root);
        Complex denominator = complex_mul_scalar(multiplicity, polynomial_complex_diff(p, multiplicity, *root));
        const Complex delta = complex_div(numerator, denominator);
        *root = complex_sub(*root, delta);
        if (complex_modulus(delta) < PRECISION) {
            break;
        }
    }
}

// Remember to free the returned vector after calling this function!
Complex_Vector polynomial_find_roots(const Vector polynomial) {
    if (polynomial.len <= 1) {
        return (Complex_Vector){0};  // Invalid operation
    }
    // Normalization of the coefficients of the polynomial
    Vector p = vector_alloc(polynomial.len);
    if (!vector_is_valid(p)) {
        return (Complex_Vector){0};
    }
    p.data[polynomial.len - 1] = 1.0;
    for (size_t i = 0; (i + 1) < polynomial.len; i++) {
        p.data[i] = polynomial.data[i] / polynomial.data[polynomial.len - 1];
    }
    Complex_Vector roots = complex_vector_alloc(polynomial.len - 1);
    if (!complex_vector_is_valid(roots)) {
        vector_dealloc(&p);
        return (Complex_Vector){0};
    }
    // Check for null roots
    size_t i = 0;
    while (are_close(p.data[0], 0.0, PRECISION)) {
        Vector division_result = (Vector){0};
        roots.data[i] = complex_init(p.data[0], 0.0);
        polynomial_ruffini_division(p, 0.0, &division_result);
        vector_dealloc(&p);
        if (!vector_is_valid(division_result)) {
            complex_vector_dealloc(&roots);
            return (Complex_Vector){0};
        }
        p = division_result;
        i++;
    }
    // Find remaining roots
    while (i < roots.len) {
        Complex root = polynomial_root_guess(p);
        uint16_t multiplicity = polynomial_find_root(p, &root);
        polynomial_root_refinement(polynomial, &root, multiplicity);
        while ((multiplicity--) && (i < roots.len)) {
            Vector division_result = (Vector){0};
            if (complex_is_real(root)) {
                polynomial_ruffini_division(p, root.real, &division_result);
                roots.data[i++] = root;
            } else {
                const double squared_modulus = root.real * root.real + root.imag * root.imag;
                Vector remainder = polynomial_quadratic_division(p, 2.0 * root.real, squared_modulus, &division_result);
                vector_dealloc(&remainder);
                roots.data[i++] = root;
                if (i > roots.len) {
                    break;
                }
                roots.data[i++] = complex_conjugate(root);
            }
            vector_dealloc(&p);
            if (!vector_is_valid(division_result)) {
                complex_vector_dealloc(&roots);
                return (Complex_Vector){0};
            }
            p = division_result;
        }
    }
    vector_dealloc(&p);
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