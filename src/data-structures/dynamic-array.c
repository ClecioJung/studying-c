#include <stdio.h>
#include <stdlib.h>

#define INITIAL_ARRAY_CAPACITY 4

typedef struct {
    int *values;
    size_t capacity;
    size_t length;
} Dynamic_Array;

void init_array(Dynamic_Array *const array) {
    array->length = 0;
    array->capacity = INITIAL_ARRAY_CAPACITY;
    array->values = malloc(array->capacity * sizeof(int));
    if (array->values == NULL) {
        array->capacity = 0;
    }
}

void free_array(Dynamic_Array *const array) {
    array->length = 0;
    array->capacity = 0;
    free(array->values);
    array->values = NULL;
}

void print_array(Dynamic_Array *const array) {
    printf("Array length: %ld, capacity: %ld\n", array->length, array->capacity);
    for (size_t i = 0; i < array->length; i++) {
        printf("[%03ld]: %d\n", i, array->values[i]);
    }
    printf("\n");
}

void grow_capacity_if_needed(Dynamic_Array *const array) {
    if (array->length + 1 > array->capacity) {
        array->capacity *= 2;
        int *new_values = realloc(array->values, array->capacity * sizeof(int));
        if (new_values != NULL) {
            array->values = new_values;
        } else {
            free_array(array);
        }
    }
}

void push(Dynamic_Array *const array, const int value) {
    grow_capacity_if_needed(array);
    array->values[array->length] = value;
    array->length++;
}

int pop(Dynamic_Array *const array) {
    array->length--;
    return array->values[array->length];
}

// Inserts element at the beginning of the array
void unshift(Dynamic_Array *const array, const int value) {
    grow_capacity_if_needed(array);
    for (size_t i = array->length; i > 0; i--) {
        array->values[i] = array->values[i - 1];
    }
    array->values[0] = value;
    array->length++;
}

// Removes and returns first element of the array
int shift(Dynamic_Array *const array) {
    if (array->length == 0) {
        return -1;
    }
    int value = array->values[0];
    array->length--;
    for (size_t i = 0; i < array->length; i++) {
        array->values[i] = array->values[i + 1];
    }
    return value;
}

int value_at(Dynamic_Array *const array, const size_t index) {
    return array->values[index];
}

size_t find_in_array(Dynamic_Array *const array, const int value) {
    // Sequential search
    for (size_t i = 0; i < array->length; i++) {
        if (array->values[i] == value) {
            return i;
        }
    }
    return array->length + 1;
}

void delete_at(Dynamic_Array *const array, const size_t index) {
    if (index >= array->length) {
        return;
    }
    array->length--;
    for (size_t i = index; i < array->length; i++) {
        array->values[i] = array->values[i + 1];
    }
}

int main(void) {
    Dynamic_Array array;
    init_array(&array);
    printf("Array before insertion:\n");
    print_array(&array);
    for (int i = 0; i < 16; i++) {
        push(&array, i);
    }
    printf("Array after insertion:\n");
    print_array(&array);

    printf("Pop a value from array: %d\n", pop(&array));
    printf("Array after pop:\n");
    print_array(&array);

    printf("Value at index 5: %d\n", value_at(&array, 5));
    printf("Find value 10 at position %ld\n", find_in_array(&array, 10));

    printf("\nArray after delete element at position 2:\n");
    delete_at(&array, 2);
    print_array(&array);

    printf("Shift a value from array: %d\n", shift(&array));
    printf("Array after shift:\n");
    print_array(&array);

    printf("Array after insertion at the beginning:\n");
    for (int i = 0; i < 5; i++) {
        unshift(&array, 100 + i);
    }
    print_array(&array);

    free_array(&array);
    return EXIT_SUCCESS;
}
