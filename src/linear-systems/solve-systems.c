#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/linear-systems.h"
#include "../../lib/scalar.h"

#define PRECISION 1e-9

typedef struct {
    System_Type type;
    Matrix A;
    Vector b;
} System_Test_Cases;

typedef System_Type (*Linear_Systems_Fn)(const Matrix A, const Vector b, Vector *const x);

typedef struct {
    const char *const name;
    const Linear_Systems_Fn method_fn;
} Linear_Systems_Algorithm;

int main(void) {
    const Linear_Systems_Algorithm algorithms[] = {
        {"Gaussian elimination", gaussian_elimination},
        {"Gauss-Jordan method", gauss_jordan},
        {"LU method", lu_solving},
    };
    System_Test_Cases system_tests[4];
    {  // System with single unique solution
        system_tests[0].type = Solvable;
        system_tests[0].A = matrix_alloc(2, 2);
        system_tests[0].b = vector_alloc(2);
        matrix_set(system_tests[0].A, 0, 0, 1.0);
        matrix_set(system_tests[0].A, 0, 1, 2.0);
        matrix_set(system_tests[0].A, 1, 0, 2.0);
        matrix_set(system_tests[0].A, 1, 1, 1.0);
        system_tests[0].b.data[0] = 3.0;
        system_tests[0].b.data[1] = 3.0;
    }
    {  // System underdetermined (not all the rows are independent)
        system_tests[1].type = Underdetermined;
        system_tests[1].A = matrix_alloc(2, 2);
        system_tests[1].b = vector_alloc(2);
        matrix_set(system_tests[1].A, 0, 0, 1.0);
        matrix_set(system_tests[1].A, 0, 1, 2.0);
        matrix_set(system_tests[1].A, 1, 0, 2.0);
        matrix_set(system_tests[1].A, 1, 1, 4.0);
        system_tests[1].b.data[0] = 3.0;
        system_tests[1].b.data[1] = 6.0;
    }
    {  // System underdetermined (more columns than rows)
        system_tests[2].type = Underdetermined;
        system_tests[2].A = matrix_alloc(2, 3);
        system_tests[2].b = vector_alloc(2);
        matrix_set(system_tests[2].A, 0, 0, 1.0);
        matrix_set(system_tests[2].A, 0, 1, 2.0);
        matrix_set(system_tests[2].A, 0, 2, 3.0);
        matrix_set(system_tests[2].A, 1, 0, 1.0);
        matrix_set(system_tests[2].A, 1, 1, 1.0);
        matrix_set(system_tests[2].A, 1, 2, 1.0);
        system_tests[2].b.data[0] = 3.0;
        system_tests[2].b.data[1] = 2.0;
    }
    {  // System overdetermined (more rows than columns)
        system_tests[3].type = Overdetermined;
        system_tests[3].A = matrix_alloc(3, 2);
        system_tests[3].b = vector_alloc(3);
        matrix_set(system_tests[3].A, 0, 0, 1.0);
        matrix_set(system_tests[3].A, 0, 1, 2.0);
        matrix_set(system_tests[3].A, 1, 0, 2.0);
        matrix_set(system_tests[3].A, 1, 1, 0.0);
        matrix_set(system_tests[3].A, 2, 0, 0.0);
        matrix_set(system_tests[3].A, 2, 1, 1.0);
        system_tests[3].b.data[0] = 3.0;
        system_tests[3].b.data[1] = 2.0;
        system_tests[3].b.data[2] = 2.0;
    }
    for (size_t i = 0; i < (sizeof(system_tests) / sizeof(system_tests[0])); i++) {
        printf("\nSolving the system nº %zu:\n", (i + 1));
        printf("Matrix A:\n");
        matrix_print(system_tests[i].A);
        printf("Vector b:\n");
        vector_print(system_tests[i].b);
        for (size_t j = 0; j < sizeof(algorithms) / sizeof(algorithms[0]); j++) {
            Vector x;
            System_Type type = algorithms[j].method_fn(system_tests[i].A, system_tests[i].b, &x);
            if (type != system_tests[i].type) {
                fprintf(stderr, "The %s couldn't get the correct type of the system!\n", algorithms[j].name);
                return EXIT_FAILURE;
            }
            if (type == Solvable) {
                if (j == 0) {
                    printf("The solution is:\n");
                    vector_print(x);
                }
                Vector result = matrix_mul_vector(system_tests[i].A, x);
                if (vector_are_equal(result, system_tests[i].b)) {
                    printf("The %s found the correct solution!\n", algorithms[j].name);
                } else {
                    fprintf(stderr, "The %s didn't found the correct solution!\n", algorithms[j].name);
                    return EXIT_FAILURE;
                }
                vector_dealloc(&result);
            } else if (type == Overdetermined) {
                if (j == 0) {
                    printf("This system is overdetermined!\n\n");
                }
            } else if (type == Underdetermined) {
                if (j == 0) {
                    printf("This system is underdetermined!\n\n");
                }
            } else {
                fprintf(stderr, "The %s produced an error!\n", algorithms[j].name);
                return EXIT_FAILURE;
            }
            vector_dealloc(&x);
        }
        matrix_dealloc(&system_tests[i].A);
        vector_dealloc(&system_tests[i].b);
    }
    return EXIT_SUCCESS;
}
