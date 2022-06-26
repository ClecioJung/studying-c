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

#ifndef __LINEAR_SYSTEMS_H
#define __LINEAR_SYSTEMS_H

#include <stdbool.h>

#include "matrix.h"
#include "vector.h"

typedef enum {
    Error = 0,        // Couldn't solve the system
    Invalid,          // System invalid (disagreeing number of rows and columns)
    Underdetermined,  // System with infinitely many solutions (Normally, a system with fewer equations than unknowns)
    Overdetermined,   // System with no solution (Normally, a system with more equations than unknowns)
    Solvable,         // Single unique solution
} System_Type;

Vector back_substitution(const Matrix A, const Vector b);
Vector forward_substitution(const Matrix A, const Vector b);
System_Type gaussian_elimination(const Matrix A, const Vector b, Vector *const x);
System_Type gauss_jordan(const Matrix A, const Vector b, Vector *const x);
System_Type lu_solving(const Matrix A, const Vector b, Vector *const x);
System_Type tridiagonal_solving(const Vector t, const Vector r, const Vector d, const Vector b, Vector *const x);
System_Type jacobi_method(const Matrix A, const Vector b, Vector *const x);
System_Type gauss_seidel(const Matrix A, const Vector b, Vector *const x);
bool columns_condition(const Matrix A);
bool rows_condition(const Matrix A);
bool sassenfeld_condition(const Matrix A);

// 'over' functions override the contents of their arguments,
// avoiding the need to allocate more memory for the results
void partial_pivoting_over(const Matrix A, const Vector b);
void back_substitution_over(const Matrix A, const Vector b);
void forward_substitution_over(const Matrix A, const Vector b);
System_Type gaussian_elimination_over(const Matrix A, const Vector b);
System_Type gauss_jordan_over(const Matrix A, const Vector b);
System_Type lu_solving_over(const Matrix A, const Vector b);
System_Type tridiagonal_solving_over(const Vector t, const Vector r, const Vector d, const Vector b);

#endif  // __LINEAR_SYSTEMS_H

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