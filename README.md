# Overview

This is a set of simple C codes developed to study algorithms and concepts of computer science. It posseses the following characteristics:

- Developed with C11;
- This code uses some UNIX functions, so it is platform dependent (it was only tested on Linux);
- In this studies, I found this references to be very usefull:
    - [Numerical Methods for Engineers](https://books.google.com.br/books?id=SA1LPgAACAAJ&hl=pt-BR&source=gbs_book_other_versions);
    - [C, the Complete Reference](https://books.google.com.br/books/about/C_the_Complete_Reference.html?id=hHc_AQAAIAAJ&redir_esc=y);
    - [Numerical Linear Algebra](https://www.google.com.br/books/edition/Numerical_Linear_Algebra/4Mou5YpRD_kC?hl=pt-BR&gbpv=1&dq=inauthor:%22Lloyd+N.+Trefethen%22&printsec=frontcover);
    - [Cálculo Numérico Computacional](https://sergiopeters.prof.ufsc.br/livro-calculo-numerico-computacional/);
- Implemented [data structures](https://en.wikipedia.org/wiki/Data_structure):
    - [Linked list](https://en.wikipedia.org/wiki/Linked_list);
    - [Doubly linked list](https://en.wikipedia.org/wiki/Doubly_linked_list);
    - [XOR linked list](https://en.wikipedia.org/wiki/XOR_linked_list);
    - [Dynamic array](https://en.wikipedia.org/wiki/Dynamic_array);
    - [Circular buffer](https://en.wikipedia.org/wiki/Circular_buffer);
    - [Binary tree](https://en.wikipedia.org/wiki/Binary_tree) (In order to study binary trees, a simple [expression parser](https://en.wikipedia.org/wiki/Parsing_expression_grammar) was developed, as an example);
- Implemented [sorting algorithms](https://en.wikipedia.org/wiki/Sorting_algorithm):
    - [Bubble sort](https://en.wikipedia.org/wiki/Bubble_sort);
    - [Select sort](https://en.wikipedia.org/wiki/Selection_sort);
    - [Insert sort](https://en.wikipedia.org/wiki/Insertion_sort);
    - [Shell sort](https://en.wikipedia.org/wiki/Shellsort);
    - [Merge sort](https://en.wikipedia.org/wiki/Merge_sort);
    - [Heap sort](https://en.wikipedia.org/wiki/Heapsort);
    - [Quicksort](https://en.wikipedia.org/wiki/Quicksort);
- Implemented [search algorithms](https://en.wikipedia.org/wiki/Search_algorithm):
    - [Sequential search](https://en.wikipedia.org/wiki/Linear_search);
    - [Binary search](https://en.wikipedia.org/wiki/Search_algorithm);
- Implemented approximations for constants:
    - [Euler's constant](https://en.wikipedia.org/wiki/Euler%27s_constant);
    - [Pi](https://en.wikipedia.org/wiki/Pi);
- Implemented scalar functions:
    - [Exponentiation by a integer positive number](https://en.wikipedia.org/wiki/Exponentiation);
    - [Square root](https://en.wikipedia.org/wiki/Square_root);
    - [N-th root](https://en.wikipedia.org/wiki/Nth_root);
    - [Exponential](https://en.wikipedia.org/wiki/Exponential_function);
    - [Trigonometric functions](https://en.wikipedia.org/wiki/Trigonometric_functions): sine, cosine, tangent and arctangent functions;
    - The [atan2](https://en.wikipedia.org/wiki/Atan2) function;
- Implemented operations with [complex numbers](https://en.wikipedia.org/wiki/Complex_number):
    - Sum andsubtraction of complex numbers;
    - Multiplication and division of complex numbers;
    - Inversion of complex numbers;
    - Modulus and argument of complex numbers;
    - Construction of complex numbers by its polar form;
    - Complex conjugate;
    - Exponentiation by a integer positive number;
- Implemented operations with [vectors](https://en.wikipedia.org/wiki/Vector_(mathematics_and_physics)):
    - Sum and subtraction of vectors;
    - Multiplication of vectors by scalars;
    - [Dot product](https://en.wikipedia.org/wiki/Dot_product);
    - [Cross product](https://en.wikipedia.org/wiki/Cross_product);
    - [Euclidean norm](https://en.wikipedia.org/wiki/Euclidean_space#Euclidean_norm);
- Implemented operations with [polynomials](https://en.wikipedia.org/wiki/Polynomial):
    - Addidtion and multiplication of polynomials;
    - Differenctiation of polynomials;
    - Evaluation with [Horner's method](https://en.wikipedia.org/wiki/Horner%27s_method);
    - Division by a binomial using [Ruffini's rule](https://en.wikipedia.org/wiki/Ruffini%27s_rule);
    - Computation of [polynomial root bounds](https://en.wikipedia.org/wiki/Geometrical_properties_of_polynomial_roots#Bounds_on_all_roots);
    - [Polynomial root finding](https://en.wikipedia.org/wiki/Root-finding_algorithms);
    - [Legendre polynomials](https://en.wikipedia.org/wiki/Legendre_polynomials);
    - [Chebyshev polynomials](https://en.wikipedia.org/wiki/Chebyshev_polynomials);
- Implemented operations with [matrices](https://en.wikipedia.org/wiki/Matrix_(mathematics)):
    - Sum and subtraction of matrices;
    - Multiplication of matrices by scalars;
    - Multiplication between matrices;
    - Multiplication between matrices and vectors;
    - [Trace](https://en.wikipedia.org/wiki/Trace_(linear_algebra));
    - [Determinant](https://en.wikipedia.org/wiki/Determinant);
    - [Transpose matrix](https://en.wikipedia.org/wiki/Transpose);
    - [Symmetric matrix](https://en.wikipedia.org/wiki/Symmetric_matrix);
    - [Skew-symmetric matrix](https://en.wikipedia.org/wiki/Skew-symmetric_matrix);
    - [LU decomposition](https://en.wikipedia.org/wiki/LU_decomposition);
    - [LU Crout decomposition](https://en.wikipedia.org/wiki/Crout_matrix_decomposition);
    - [QR decomposition](https://en.wikipedia.org/wiki/QR_decomposition);
    - [Householder matrix computation](https://en.wikipedia.org/wiki/Householder_transformation#Householder_matrix);
    - [Hessenberg matrix computation](https://en.wikipedia.org/wiki/Hessenberg_matrix);
    - [Schur decomposition](https://en.wikipedia.org/wiki/Schur_decomposition) using the [QR algorithm](https://en.wikipedia.org/wiki/QR_algorithm);
    - [Eigenvalues computation](https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors);
    - [Cholesky decomposition](https://en.wikipedia.org/wiki/Cholesky_decomposition);
    - Greatest eigenvalue computation (and its eigenvector) using the [power method](https://en.wikipedia.org/wiki/Power_iteration);
    - Lowest eigenvalue computation (and its eigenvector) using the [inverse power method](https://en.wikipedia.org/wiki/Inverse_iteration);
    - [Inverse matrix](https://en.wikipedia.org/wiki/Invertible_matrix);
    - [Pseudo-inverse matrix](https://en.wikipedia.org/wiki/Generalized_inverse);
    - [Jacobian matrix](https://en.wikipedia.org/wiki/Jacobian_matrix_and_determinant);
    - [Vandermonde matrix](https://en.wikipedia.org/wiki/Vandermonde_matrix);
- Implemented algorithm for solving [systems of linear equations](https://en.wikipedia.org/wiki/System_of_linear_equations):
    - [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination);
    - Gauss-Jordan elimination;
    - [Solving linear equations using LU decomposition](https://en.wikipedia.org/wiki/LU_decomposition#Solving_linear_equations);
    - [Jacobi method](https://en.wikipedia.org/wiki/Jacobi_method);
    - [Gauss-Seidel method](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method);
    - Algorithm for solving [tridiagonal systems](https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm);
- Implemented algorithms for [finding roots of equations](https://en.wikipedia.org/wiki/Root-finding_algorithms):
    - [Bisection method](https://en.wikipedia.org/wiki/Bisection_method);
    - [Fake-position method](https://en.wikipedia.org/wiki/Regula_falsi);
    - [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method);
    - [Secant method](https://en.wikipedia.org/wiki/Secant_method);
    - [Muller's method](https://en.wikipedia.org/wiki/Muller%27s_method);
    - Newton-Raphson method for multivariable systems;
    - [Broyden's method for multivariable systems](https://en.wikipedia.org/wiki/Broyden%27s_method);
- Implemented algorithms for [curve-fitting](https://en.wikipedia.org/wiki/Curve_fitting):
    - [Polynomial interpolation](https://en.wikipedia.org/wiki/Polynomial_interpolation);
    - [Lagrange interpolation](https://en.wikipedia.org/wiki/Lagrange_polynomial);
    - [Gregory-Newton intepolation](https://en.wikipedia.org/wiki/Newton_polynomial);
    - [Linear regression](https://en.wikipedia.org/wiki/Linear_regression);
    - [Polynomail regression](https://en.wikipedia.org/wiki/Polynomial_regression);
    - Computation of [Coefficient of determination](https://en.wikipedia.org/wiki/Coefficient_of_determination) for the polynomial regression;
    - [Cubic spline interpolation](https://en.wikipedia.org/wiki/Spline_interpolation);
- Implemented algorithms for [numerical integration](https://en.wikipedia.org/wiki/Numerical_integration):
    - [Trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule);
    - [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule);
    - [Simpson's 3/8 rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_3/8_rule);
    - [Gauss-Legendre quadrature](https://en.wikipedia.org/wiki/Gauss%E2%80%93Legendre_quadrature);
- Implemented algorithms for solving [Ordinary Differential Equations (ODE)](https://en.wikipedia.org/wiki/Ordinary_differential_equation):
    - [Euler method](https://en.wikipedia.org/wiki/Euler_method);
    - [Heun's method](https://en.wikipedia.org/wiki/Heun%27s_method);
    - Third order Runge-Kutta method;
    - [Fourth order Runge-Kutta method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods);
    - Butcher’s (1964) fifth-order RK method;

# Usage

Download this project and compile it by typing the command `make` in its folder. Next, just run one of the executables located in the `bin` folder. Here is an example:

```console
$ make
$ ./bin/binary-tree
```

In order to check if there is any memory leak, use the `valgrind` command:

```console
$ valgrind --tool=memcheck ./bin/binary-tree
```

To run all the tests, use the following command:

```console
$ make test
```