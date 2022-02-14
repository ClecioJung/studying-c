#include <stdio.h>
#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    const size_t len = 30;
    Vector vec = random_vector(len, 0.0, (data_type)len);
    print_vector(vec);
    const data_type value = 15.0;
    const size_t index = sequential_search(vec, value);
    if (index < len) {
        printf("Sequential search: found value %g at index %ld\n", vec.data[index], index);
    } else {
        printf("Sequential search: didn't found value %g\n", value);
    }
    free_vector(&vec);
    return EXIT_SUCCESS;
}