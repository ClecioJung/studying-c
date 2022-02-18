#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define PRECISION 1e-10

typedef struct {
    const double input;
    const double expected;
} Test_Values;

int main(void) {
    const Test_Values values[] = {
        {-9.0, 0.00012340980408667956},
        {-8.0, 0.00033546262790251185},
        {-7.0, 0.0009118819655545162},
        {-6.0, 0.0024787521766663585},
        {-5.0, 0.006737946999085467},
        {-4.0, 0.01831563888873418},
        {-3.0, 0.049787068367863944},
        {-2.0, 0.1353352832366127},
        {-1.0, 0.36787944117144233},
        {0.0, 1},
        {1.0, 2.718281828459045},
        {2.0, 7.38905609893065},
        {3.0, 20.085536923187668},
        {4.0, 54.598150033144236},
        {5.0, 148.4131591025766},
        {6.0, 403.4287934927351},
        {7.0, 1096.6331584284585},
        {8.0, 2980.9579870417283},
        {9.0, 8103.083927575384},
    };
    printf("Exponential:\n");
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        const double result = exponential(values[i].input);
        if (are_close(result, values[i].expected, PRECISION)) {
            printf("exponential(%lg) = %lg\n", values[i].input, result);
        } else {
            fprintf(stderr, "Error: imprecise result. exponential(%lg) = %lg\n", values[i].input, result);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}