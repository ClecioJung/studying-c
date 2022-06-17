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

#include "complex.h"

#include <math.h>

#include "scalar.h"

#define MAX_ITERATIONS 10000
#define PRECISION 1e-10

Complex complex_init(const double real, const double imag) {
    return (Complex){
        .real = real,
        .imag = imag,
    };
}

Complex complex_polar(const double modulus, const double phase) {
    const double real = modulus * cosine(phase);
    const double imag = modulus * sine(phase);
    return complex_init(real, imag);
}

double complex_modulus(const Complex c) {
    return square_root(c.real * c.real + c.imag * c.imag);
}

double complex_argument(const Complex c) {
    return atan2(c.imag, c.real);
}

Complex complex_conjugate(const Complex c) {
    return complex_init(c.real, -c.imag);
}

Complex complex_inverse(const Complex c) {
    const double squared_modulus = c.real * c.real + c.imag * c.imag;
    const double real = c.real / squared_modulus;
    const double imag = -c.imag / squared_modulus;
    return complex_init(real, imag);
}

Complex complex_sum(const Complex a, const Complex b) {
    return complex_init(a.real + b.real, a.imag + b.imag);
}

Complex complex_sub(const Complex a, const Complex b) {
    return complex_init(a.real - b.real, a.imag - b.imag);
}

Complex complex_mul(const Complex a, const Complex b) {
    const double real = a.real * b.real - a.imag * b.imag;
    const double imag = a.real * b.imag + a.imag * b.real;
    return complex_init(real, imag);
}

Complex complex_div(const Complex a, const Complex b) {
    return complex_mul(a, complex_inverse(b));
}

bool complex_are_equal(const Complex a, const Complex b) {
    return (are_close(a.real, b.real, PRECISION) && are_close(a.imag, b.imag, PRECISION));
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