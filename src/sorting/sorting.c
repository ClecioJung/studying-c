#include "../../lib/sorting.h"

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

typedef void (*Sorting_Fn)(const Vector);

typedef struct {
    const char *const name;
    const Sorting_Fn sorting_fn;
} Sorting_Algorithm;

int test_algorithm(const Vector vec, const Sorting_Algorithm algorithm) {
    Vector copied_vec = vector_copy(vec);
    struct timeval stop, start;
    gettimeofday(&start, NULL);
    algorithm.sorting_fn(copied_vec);
    gettimeofday(&stop, NULL);
    const bool is_sorted = vector_is_sorted(copied_vec);
    if (is_sorted) {
        float delta_ms = ((float)((stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec))) / 1000.0;
        printf("Took %06.2f ms for the %s algorithm to sort the provided vector!\n", delta_ms, algorithm.name);
    } else {
        fprintf(stderr, "The algorithm %s couldn't sort the provided vector!\n", algorithm.name);
    }
    vector_dealloc(&copied_vec);
    return is_sorted ? EXIT_SUCCESS : EXIT_FAILURE;
}

int main(void) {
    const Sorting_Algorithm algorithms[] = {
        {"buble-sort", bubble_sort},
        {"select-sort", select_sort},
        {"insert-sort", insert_sort},
        {"shell-sort", shell_sort},
        {"merge-sort", merge_sort},
        {"heap-sort", heap_sort},
        {"quicksort", quicksort},
    };
    const size_t len = 10000;
    Vector vec = vector_random(len, 0.0, 10.0 * len);
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(algorithms[0]); i++) {
        if (test_algorithm(vec, algorithms[i]) != EXIT_SUCCESS) {
            return EXIT_FAILURE;
        }
    }
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}
