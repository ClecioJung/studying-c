#include <stdlib.h>

#define MATRIX_IMPLEMENTATION
#include "../math-vector.h"

int main(void) {
    const size_t len = 30;
    Vector vec = random_vector(len, 0.0, 10.0 * len);
    bubble_sort(vec);
    print_vector(vec);
    free_vector(&vec);
    return EXIT_SUCCESS;
}
