#include <stdio.h>
#include <stdlib.h>

#include "../../lib/curve-fitting.h"
#include "../../lib/scalar.h"

#define PRECISION 1e-12

int main(void) {
    Vector x = vector_alloc(5);
    Vector y = vector_alloc(5);
    x.data[0] = 0.0;
    x.data[1] = 1.0;
    x.data[2] = 2.0;
    x.data[3] = 3.0;
    x.data[4] = 4.0;
    y.data[0] = -3.0;
    y.data[1] = -2.0;
    y.data[2] = 5.0;
    y.data[3] = 24.0;
    y.data[4] = 61.0;
    const double value = 2.5;
    printf("Considering the data:\nx =\n");
    vector_print(x);
    printf("y =\n");
    vector_print(y);
    {
        printf("The natural spline interpolation resulted in:\n");
        Spline spline = spline_interpolation(x, y, Spline_Natural, 0.0, 0.0);
        spline_print(spline);
        for (size_t i = 0; i < x.len; i++) {
            const double result = spline_evaluation(spline, x.data[i]);
            if (!are_close(result, y.data[i], PRECISION)) {
                fprintf(stderr, "Spline interpolation: couldn't interpolate correctly the provided values!\n");
                return EXIT_FAILURE;
            }
        }
        printf("The spline interpolation avaliated at %lg results in %lg\n\n", value, spline_evaluation(spline, value));
        spline_dealloc(&spline);
    }
    {
        printf("The quadratic spline interpolation resulted in:\n");
        Spline spline = spline_interpolation(x, y, Spline_Quadratic, 0.0, 0.0);
        spline_print(spline);
        for (size_t i = 0; i < x.len; i++) {
            const double result = spline_evaluation(spline, x.data[i]);
            if (!are_close(result, y.data[i], PRECISION)) {
                fprintf(stderr, "Spline interpolation: couldn't interpolate correctly the provided values!\n");
                return EXIT_FAILURE;
            }
        }
        printf("The spline interpolation avaliated at %lg results in %lg\n\n", value, spline_evaluation(spline, value));
        spline_dealloc(&spline);
    }
    {
        printf("The clamped spline interpolation resulted in:\n");
        Spline spline = spline_interpolation(x, y, Spline_Clamped, 0.0, 24.0);
        spline_print(spline);
        for (size_t i = 0; i < x.len; i++) {
            const double result = spline_evaluation(spline, x.data[i]);
            if (!are_close(result, y.data[i], PRECISION)) {
                fprintf(stderr, "Spline interpolation: couldn't interpolate correctly the provided values!\n");
                return EXIT_FAILURE;
            }
        }
        printf("The spline interpolation avaliated at %lg results in %lg\n\n", value, spline_evaluation(spline, value));
        spline_dealloc(&spline);
    }
    vector_dealloc(&x);
    vector_dealloc(&y);
    return EXIT_SUCCESS;
}