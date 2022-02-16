# Overview

This is a set of simple C codes developed to study algorithms and concepts of computer science. It posseses the following characteristics:

- Developed with C11;
- This code uses some UNIX functions, so it is platform dependent (it was tested only on Linux);
- In this studies, I found this references to be very usefull:
    - [Numerical Methods for Engineers](https://books.google.com.br/books?id=SA1LPgAACAAJ&hl=pt-BR&source=gbs_book_other_versions);
    - [C, the Complete Reference](https://books.google.com.br/books/about/C_the_Complete_Reference.html?id=hHc_AQAAIAAJ&redir_esc=y);
- In the development of this code I created some very useful [stb-style header-only librarys](https://github.com/nothings/stb), such as:
    - [math-vector.h](./src/math-vector.h): Useful math functions for dealing with scalars, vectors, matrices, solve linear systems, curve fitting and more. Also includes functions for sorting and searchnig vectors;
    - [roots.h](./src/roots-of-eqs/roots.h): Useful functions for finding roots of equations;
    - [integration.h](./src/integration/integration.h): Useful functions for numerical integration;
    - [ode.h](./src/diff-eqs/ode.h): Useful functions for solving Ordinary Differential Equations (ODE) numerically;
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
    - [Square root](https://en.wikipedia.org/wiki/Square_root);
    - [Exponential](https://en.wikipedia.org/wiki/Exponential_function);
- Implemented operations with [vectors](https://en.wikipedia.org/wiki/Vector_(mathematics_and_physics)):
    - Sum/subtraction of vectors and multiplication by scalars;
    - [Dot product](https://en.wikipedia.org/wiki/Dot_product);
    - [Cross product](https://en.wikipedia.org/wiki/Cross_product);
    - [Euclidean norm](https://en.wikipedia.org/wiki/Euclidean_space#Euclidean_norm);
- Implemented operations with [matrices](https://en.wikipedia.org/wiki/Matrix_(mathematics)):
    - Sum/subtraction of matrices and multiplication by scalars;
    - Multiplication between matrices, and matrices by vectors;
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
    - Greatest eigenvalue computation (and its eigenvector) using the [power method](https://en.wikipedia.org/wiki/Power_iteration);
    - [Inverse matrix](https://en.wikipedia.org/wiki/Invertible_matrix);
    - [Pseudo-inverse matrix](https://en.wikipedia.org/wiki/Generalized_inverse);
- Implemented algorithm for solving [systems of linear equations](https://en.wikipedia.org/wiki/System_of_linear_equations):
    - [Gaussian elimination](https://en.wikipedia.org/wiki/Gaussian_elimination);
    - Gauss-Jordan elimination;
    - [Solving linear equations using LU decomposition](https://en.wikipedia.org/wiki/LU_decomposition#Solving_linear_equations);
    - [Jacobi method](https://en.wikipedia.org/wiki/Jacobi_method);
    - [Gauss-Seidel method](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method);
- Implemented algorithms for [finding roots of equations](https://en.wikipedia.org/wiki/Root-finding_algorithms):
    - [Bisection method](https://en.wikipedia.org/wiki/Bisection_method);
    - [Fake-position method](https://en.wikipedia.org/wiki/Regula_falsi);
    - [Newton-Raphson method](https://en.wikipedia.org/wiki/Newton%27s_method);
    - [Secant method](https://en.wikipedia.org/wiki/Secant_method);
- Implemented algorithms for [curve-fitting](https://en.wikipedia.org/wiki/Curve_fitting):
    - [Lagrange interpolation](https://en.wikipedia.org/wiki/Lagrange_polynomial);
    - [Linear regression](https://en.wikipedia.org/wiki/Linear_regression);
    - [Polynomail regression](https://en.wikipedia.org/wiki/Polynomial_regression);
- Implemented algorithms for [numerical integration](https://en.wikipedia.org/wiki/Numerical_integration):
    - [Trapezoidal rule](https://en.wikipedia.org/wiki/Trapezoidal_rule);
    - [Simpson's rule](https://en.wikipedia.org/wiki/Simpson%27s_rule);
    - [Simpson's 3/8 rule](https://en.wikipedia.org/wiki/Simpson%27s_rule#Simpson's_3/8_rule);
- Implemented algorithms for solving [Ordinary Differential Equations (ODE)](https://en.wikipedia.org/wiki/Ordinary_differential_equation):
    - [Euler method](https://en.wikipedia.org/wiki/Euler_method);
    - [Heun's method](https://en.wikipedia.org/wiki/Heun%27s_method);
    - Third order Runge-Kutta method;
    - [Fourth order Runge-Kutta method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta_methods);
    - Butcher’s (1964) fifth-order RK method;
# Usage

Download this project and compile it by typing the command `make` in its folder. Next, just run one of the executables. Here is an example:

```console
$ make
$ ./linked-list
$ ./doubly-linked-list
$ ./xor-linked-list
$ ./dynamic-array
$ ./circular-buffer
$ ./binary-tree
$ ./bubble-sort
$ ./select-sort
$ ./insert-sort
$ ./shell-sort
$ ./merge-sort
$ ./heap-sort
$ ./quicksort
$ ./compare-sorting
$ ./sequential-search
$ ./binary-search
$ ./compare-search
$ ./euler
$ ./pi
$ ./dot-product
$ ./cross-product
$ ./euclidean-norm
$ ./matrix-transpose
$ ./determinant
$ ./lu-decomposition
$ ./lu-crout-decomposition
$ ./matrix-inverse
$ ./pseudo-inverse
$ ./gaussian-elimination
$ ./gauss-jordan
$ ./lu-solving
$ ./jacobi-method
$ ./gauss-seidel
$ ./roots
$ ./lagrange-interpolation
$ ./linear-regression
$ ./polynomial-regression
$ ./integration
$ ./ode
$ ./square-root
$ ./exponential
$ ./trace
$ ./matrix-symmetric
$ ./qr-decomposition
$ ./householder-matrix
$ ./upper-hessenberg-matrix
$ ./schur-decomposition
$ ./eigenvalues
$ ./power-method
```

In order to check if there is any memory leak, run the following commands (it is necessary to have `valgrind` installed):

```console
$ valgrind --tool=memcheck ./linked-list
$ valgrind --tool=memcheck ./doubly-linked-list
$ valgrind --tool=memcheck ./xor-linked-list
$ valgrind --tool=memcheck ./dynamic-array
$ valgrind --tool=memcheck ./binary-tree
```