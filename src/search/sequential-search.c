#include <stdio.h>
#include <stdlib.h>

#include "../../lib/search.h"

int main(void) {
    const size_t len = 30;
    Vector vec = vector_random(len, 0.0, (double)len);
    vector_print(vec);
    const double value = 15.0;
    const size_t index = sequential_search(vec, value);
    if (index < len) {
        printf("Sequential search: found value %g at index %zu\n", vec.data[index], index);
    } else {
        printf("Sequential search: didn't found value %g\n", value);
    }
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}