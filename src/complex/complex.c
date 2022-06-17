#include "../../lib/complex.h"

#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define PRECISION 1e-10

int main(void) {
    Complex a = (Complex){.real = 1.0, .imag = 0.0};
    Complex b = complex_init(square_root(2.0), square_root(2.0));
    Complex c = complex_polar(1.0, 0.0);
    Complex b_inv = complex_inverse(b);
    {
        if (are_close(complex_modulus(a), 1.0, PRECISION)) {
            printf("The modulus of the complex number 'a' is equal to 1.0!\n");
        } else {
            fprintf(stderr, "The modulus of the complex number 'a' should've been equal to 1.0!\n");
            return EXIT_FAILURE;
        }
        if (are_close(complex_modulus(b), 2.0, PRECISION)) {
            printf("The modulus of the complex number 'b' is equal to 2.0!\n");
        } else {
            fprintf(stderr, "The modulus of the complex number 'b' should've been equal to 2.0!\n");
            return EXIT_FAILURE;
        }
        if (are_close(complex_modulus(c), 1.0, PRECISION)) {
            printf("The modulus of the complex number 'c' is equal to 1.0!\n");
        } else {
            fprintf(stderr, "The modulus of the complex number 'c' should've been equal to 1.0!\n");
            return EXIT_FAILURE;
        }
    }
    {
        if (are_close(complex_argument(a), 0.0, PRECISION)) {
            printf("The argument of the complex number 'a' is equal to 0.0!\n");
        } else {
            fprintf(stderr, "The argument of the complex number 'a' should've been equal to 0.0!\n");
            return EXIT_FAILURE;
        }
        if (are_close(complex_argument(b), PI / 4.0, PRECISION)) {
            printf("The argument of the complex number 'b' is equal to pi/2!\n");
        } else {
            fprintf(stderr, "The argument of the complex number 'b' should've been equal to pi/4!\n");
            return EXIT_FAILURE;
        }
        if (are_close(complex_argument(c), 0.0, PRECISION)) {
            printf("The argument of the complex number 'c' is equal to 0.0!\n");
        } else {
            fprintf(stderr, "The argument of the complex number 'c' should've been equal to 0.0!\n");
            return EXIT_FAILURE;
        }
    }
    {
        if (complex_are_equal(b, complex_polar(complex_modulus(b), complex_argument(b)))) {
            printf("Could correctly recreate the complex number 'b' from its polar form!\n");
        } else {
            fprintf(stderr, "Couldn't correctly recreate the complex number 'b' from its polar form!\n");
            return EXIT_FAILURE;
        }
    }
    {
        if (complex_are_equal(a, c)) {
            printf("The complex numbers 'a' and 'c' are equal!\n");
        } else {
            fprintf(stderr, "The complex numbers 'a' and 'c' should've been equal!\n");
            return EXIT_FAILURE;
        }
    }
    {
        if (complex_are_equal(a, complex_mul(b, b_inv)) && are_close(complex_modulus(b_inv), 0.5, PRECISION)) {
            printf("The inverse of 'b' was correctly calculated!\n");
        } else {
            fprintf(stderr, "The inverse of 'b' was NOT correctly calculated!\n");
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}