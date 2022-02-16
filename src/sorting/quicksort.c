#include <stdlib.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    const size_t len = 30;
    Vector vec = vector_random(len, 0.0, 10.0 * len);
    quicksort(vec);
    vector_print(vec);
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}
