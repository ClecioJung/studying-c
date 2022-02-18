#include <stdio.h>
#include <stdlib.h>

#include "../../lib/search.h"
#include "../../lib/sorting.h"

int main(void) {
    const size_t len = 30;
    Vector vec = vector_random(len, 0.0, (double)len);
    quicksort(vec);
    vector_print(vec);
    const double value = 15.0;
    const size_t index = binary_search(vec, value);
    if (index < len) {
        printf("Binary search: found value %g at index %zu\n", vec.data[index], index);
    } else {
        printf("Binary search: didn't found value %g\n", value);
    }
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}
