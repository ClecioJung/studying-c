#include <stdio.h>
#include <stdlib.h>

#include "../../lib/scalar.h"
#include "../../lib/search.h"

#define PRECISION 1e-10

int main(void) {
    const size_t len = 30;
    Vector vec = vector_alloc(len);
    for (size_t i = 0; i < vec.len; i++) {
        vec.data[i] = (double)i;
    }
    vector_print(vec);
    const double value = 15.0;
    const size_t index = sequential_search(vec, value);
    if ((index < len) && are_close(value, vec.data[index], PRECISION)) {
        printf("Sequential search: found value %lg at index %zu\n", vec.data[index], index);
    } else {
        fprintf(stderr, "Sequential search: didn't found value %lg\n", value);
        return EXIT_FAILURE;
    }
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}