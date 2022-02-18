#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#include "../../lib/scalar.h"
#include "../../lib/search.h"
#include "../../lib/sorting.h"

#define PRECISION 1e-10

typedef size_t (*Search_Fn)(const Vector, const double);

typedef struct {
    const char *const name;
    const Search_Fn search_fn;
} Search_Algorithm;

int test_algorithm(const Vector vec, const Search_Algorithm algorithm) {
    struct timeval stop, start;
    const double value = ((double)vec.len) / 2.0;
    gettimeofday(&start, NULL);
    algorithm.search_fn(vec, value);
    gettimeofday(&stop, NULL);
    unsigned long int delta_us = (stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec);
    printf("Took %04ld us for the %s algorithm to search the provided vector!\n", delta_us, algorithm.name);
    return EXIT_SUCCESS;
}

int main(void) {
    const Search_Algorithm algorithms[] = {
        {"sequential-search", sequential_search},
        {"binary-search", binary_search},
    };
    const size_t len = 1000000;
    Vector vec = vector_random(len, 0.0, (double)len);
    quicksort(vec);
    for (size_t i = 0; i < sizeof(algorithms) / sizeof(algorithms[0]); i++) {
        test_algorithm(vec, algorithms[i]);
    }
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}
