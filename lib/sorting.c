// MIT License

// Copyright (c) 2022 CLECIO JUNG <clecio.jung@gmail.com>

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

//------------------------------------------------------------------------------
// SOURCE
//------------------------------------------------------------------------------

#include "sorting.h"

#include <stdbool.h>

#include "scalar.h"

bool vector_is_sorted(const Vector vec) {
    for (size_t i = 0; (i + 1) < vec.len; i++) {
        if (vec.data[i] > vec.data[i + 1]) {
            return false;
        }
    }
    return true;
}

void bubble_sort(const Vector vec) {
    for (size_t i = 1; i < vec.len; i++) {
        for (size_t j = (vec.len - 1); j >= i; j--) {
            if (vec.data[j - 1] > vec.data[j]) {
                swap(&vec.data[j - 1], &vec.data[j]);
            }
        }
    }
}

void select_sort(const Vector vec) {
    for (size_t i = 0; (i + 1) < vec.len; i++) {
        size_t k = i;
        double temp = vec.data[i];
        for (size_t j = (i + 1); j < vec.len; j++) {
            if (vec.data[j] < temp) {
                k = j;
                temp = vec.data[j];
            }
        }
        if (k != i) {
            vec.data[k] = vec.data[i];
            vec.data[i] = temp;
        }
    }
}

void insert_sort(const Vector vec) {
    for (size_t i = 1; i < vec.len; i++) {
        const double temp = vec.data[i];
        size_t j = (i - 1);
        while (j < vec.len && temp < vec.data[j]) {
            vec.data[j + 1] = vec.data[j];
            j--;
        }
        vec.data[j + 1] = temp;
    }
}

void shell_sort(const Vector vec) {
    // Ciura gap sequence
    static const size_t gaps[] = {701, 301, 132, 57, 23, 10, 4, 1};
    static const size_t gap_len = sizeof(gaps) / sizeof(gaps[0]);

    for (size_t k = 0; k < gap_len; k++) {
        const size_t gap = gaps[k];
        for (size_t i = gap; i < vec.len; i++) {
            const double temp = vec.data[i];
            size_t j = (i - gap);
            while (j < vec.len && temp < vec.data[j]) {
                vec.data[j + gap] = vec.data[j];
                j = (j - gap);
            }
            vec.data[j + gap] = temp;
        }
    }
}

// Heapify a subtree rooted with node in array at index
static void heapify(double *const array, const size_t len, const size_t index) {
    size_t largest = index;
    const size_t left = 2 * index + 1;
    const size_t right = 2 * index + 2;
    // If left child is larger than root
    if (left < len && array[left] > array[largest]) {
        largest = left;
    }
    // If right child is larger than largest so far
    if (right < len && array[right] > array[largest]) {
        largest = right;
    }
    // If largest is not root
    if (largest != index) {
        swap(&array[largest], &array[index]);
        // Recursively heapify the affected sub-tree
        heapify(array, len, largest);
    }
}

void heap_sort(const Vector vec) {
    // Build heap (rearrange array)
    for (size_t i = vec.len / 2 - 1; i < vec.len; i--) {
        heapify(vec.data, vec.len, i);
    }
    // One by one extract an element from heap
    for (size_t i = vec.len - 1; i > 0; i--) {
        swap(&vec.data[0], &vec.data[i]);
        // Call max heapify on the reduced heap
        heapify(vec.data, i, 0);
    }
}

// Merge two subarrays L and M into A
static void merge(double *const array, const size_t p, const size_t q, const size_t r) {
    // Create L <- A[p..q] and M <- A[q+1..r]
    const size_t n1 = q - p + 1;
    const size_t n2 = r - q;
    double L[n1], M[n2];
    for (size_t i = 0; i < n1; i++) {
        L[i] = array[p + i];
    }
    for (size_t j = 0; j < n2; j++) {
        M[j] = array[q + 1 + j];
    }
    // Maintain current index of sub-arrays and main array
    size_t i = 0;
    size_t j = 0;
    size_t k = p;
    // Until we reach either end of either L or M, pick larger among
    // elements L and M and place them in the correct position at A[p..r]
    while (i < n1 && j < n2) {
        if (L[i] <= M[j]) {
            array[k] = L[i];
            i++;
        } else {
            array[k] = M[j];
            j++;
        }
        k++;
    }
    // When we run out of elements in either L or M,
    // pick up the remaining elements and put in A[p..r]
    while (i < n1) {
        array[k] = L[i];
        i++;
        k++;
    }
    while (j < n2) {
        array[k] = M[j];
        j++;
        k++;
    }
}

// Divide the array into two subarrays, sort them and merge them
static void split_and_sort(double *const array, const size_t left, const size_t right) {
    if (left < right) {
        double middle = left + (right - left) / 2;
        split_and_sort(array, left, middle);
        split_and_sort(array, middle + 1, right);
        merge(array, left, middle, right);
    }
}

void merge_sort(const Vector vec) {
    if (vec.len > 0) {
        split_and_sort(vec.data, 0, vec.len - 1);
    }
}

static void qs(double *const array, const size_t left, const size_t right) {
    const double middle = array[(left + right) / 2];
    size_t i = left;
    size_t j = right;
    do {
        while (array[i] < middle && i < right) {
            i++;
        }
        while (middle < array[j] && j > left && j != 0) {
            j--;
        }
        if (i <= j) {
            swap(&array[i], &array[j]);
            i++;
            if (j != 0) {
                j--;
            }
        }
    } while (i <= j);
    if (left < j) {
        qs(array, left, j);
    }
    if (i < right) {
        qs(array, i, right);
    }
}

void quicksort(const Vector vec) {
    if (vec.len > 0) {
        qs(vec.data, 0, vec.len - 1);
    }
}

//------------------------------------------------------------------------------
// END
//------------------------------------------------------------------------------

// MIT License

// Copyright (c) 2022 CLECIO JUNG <clecio.jung@gmail.com>

// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:

// The above copyright notice and this permission notice shall be included in all
// copies or substantial portions of the Software.

// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.