#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>

#define MATH_VECTOR_IMPLEMENTATION
#include "../math-vector.h"

typedef size_t (*Search_Fn)(const Vector, const data_type);

void test_algorithm(const Vector vec, const char *const algorithm, const Search_Fn search_fn) {
    struct timeval stop, start;
    const data_type value = ((data_type)vec.len) / 2.0;
    gettimeofday(&start, NULL);
    search_fn(vec, value);
    gettimeofday(&stop, NULL);
    unsigned long int delta_us = (stop.tv_sec - start.tv_sec) * 1000000 + (stop.tv_usec - start.tv_usec);
    printf("Took %04ld us for the %s algorithm to search the provided vector!\n", delta_us, algorithm);
}

int main(void) {
    const size_t len = 10000000;
    Vector vec = vector_random(len, 0.0, (data_type)len);
    quicksort(vec);
    test_algorithm(vec, "sequential-search", sequential_search);
    test_algorithm(vec, "binary-search", binary_search);
    vector_dealloc(&vec);
    return EXIT_SUCCESS;
}
