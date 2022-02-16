#include <stdlib.h>
#include <sys/time.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

typedef void (*Sorting_Fn)(const Vector);

void test_algorithm(const Vector vec, const char *const algorithm, const Sorting_Fn sorting_fn) {
    Vector copied_vec = vector_copy(vec);
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    sorting_fn(copied_vec);
    gettimeofday(&stop, NULL);
    if (vector_is_sorted(copied_vec)) {
        float delta_ms = ((float)((stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec))) / 1000.0;
        printf("Took %06.2f ms for the %s algorithm to sort the provided vector!\n", delta_ms, algorithm);
    } else {
        printf("%s: Couldn't sort the provided vector!\n", algorithm);
    }
    vector_dealloc(&copied_vec);
}

int main(void) {
    const size_t len = 100000;
    Vector vec = vector_random(len, 0.0, 10.0 * len);
    test_algorithm(vec, "buble-sort", bubble_sort);
    test_algorithm(vec, "select-sort", select_sort);
    test_algorithm(vec, "insert-sort", insert_sort);
    test_algorithm(vec, "shell-sort", shell_sort);
    test_algorithm(vec, "merge-sort", merge_sort);
    test_algorithm(vec, "heap-sort", heap_sort);
    test_algorithm(vec, "quicksort", quicksort);
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}
