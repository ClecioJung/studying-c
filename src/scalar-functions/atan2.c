#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"

#define PRECISION 1e-10

typedef struct {
    const double y;
    const double x;
    const double expected;
} Test_Values;

int main(void) {
    const Test_Values values[] = {
        {-1.0, -1.0, -2.356194490192345},
        {-1.0, 1.0, -0.7853981633974483},
        {1.0, -1.0, 2.356194490192345},
        {1.0, 1.0, 0.7853981633974483},
        {-0.5, -0.5, -2.356194490192345},
        {-0.5, 0.5, -0.7853981633974483},
        {0.5, -0.5, 2.356194490192345},
        {0.5, 0.5, 0.7853981633974483},
        {-1.0, -0.5, -2.0344439357957027},
        {-1.0, 0.5, -1.1071487177940904},
        {1.0, -0.5, 2.0344439357957027},
        {1.0, 0.5, 1.1071487177940904},
        {0.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, -1.0, 3.141592653589793},
        {1.0, 0.0, 1.5707963267948966},
        {-1.0, 0.0, -1.5707963267948966},
    };
    printf("atan2:\n");
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        const double result = atan2(values[i].y, values[i].x);
        if (are_close(result, values[i].expected, PRECISION)) {
            printf("atan2(%lg, %lg) = %lg\n", values[i].y, values[i].x, result);
        } else {
            fprintf(stderr, "Error: imprecise result. atan2(%lg, %lg) = %lg\n", values[i].y, values[i].x, result);
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}