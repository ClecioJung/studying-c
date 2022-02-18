#include <stdio.h>
#include <stdlib.h>

#include "../../lib/sorting.h"

int main(void) {
    const size_t len = 30;
    Vector vec = vector_random(len, 0.0, 10.0 * len);
    select_sort(vec);
    vector_print(vec);
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}
