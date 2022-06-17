#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "../../lib/complex.h"

typedef struct {
    const Complex base;
    const uint64_t expoent;
    const Complex expected;
} Test_Values;

int main(void) {
    const Test_Values values[] = {
        {complex_init(0.0, 0.0), 0, complex_init(1.0, 0.0)},
        {complex_init(10.0, 10.0), 0, complex_init(1.0, 0.0)},
        {complex_init(0.0, 0.0), 10, complex_init(0.0, 0.0)},
        {complex_init(-1.0, 1.0), 1, complex_init(-1.0, 1.0)},
        {complex_init(-1.0, 1.0), 2, complex_init(0.0, -2.0)},
        {complex_init(0.0, 1.0), 2, complex_init(-1.0, 0.0)},
        {complex_init(0.5, 0.5), 2, complex_init(0.0, 0.5)},
        {complex_init(0.0, 2.0), 2, complex_init(-4.0, 0.0)},
        {complex_init(1.0, 1.0), 2, complex_init(0.0, 2.0)},
        {complex_init(1.0, 1.0), 3, complex_init(-2.0, 2.0)},
        {complex_init(1.0, 1.0), 4, complex_init(-4.0, 0.0)},
        {complex_init(1.0, 1.0), 5, complex_init(-4.0, -4.0)},
        {complex_init(1.0, 1.0), 10, complex_init(0.0, 32.0)},
    };
    printf("Power:\n");
    for (size_t i = 0; i < sizeof(values) / sizeof(values[0]); i++) {
        const Complex result = complex_power(values[i].base, values[i].expoent);
        if (complex_are_equal(result, values[i].expected)) {
            printf("power(");
            complex_print(values[i].base);
            printf(", %ld) = ", values[i].expoent);
            complex_print(result);
            printf("\n");
        } else {
            fprintf(stderr, "Error: imprecise result.\n");
            printf("power(");
            complex_print(values[i].base);
            printf(", %ld) = ", values[i].expoent);
            complex_print(result);
            printf("\n");
            return EXIT_FAILURE;
        }
    }
    return EXIT_SUCCESS;
}