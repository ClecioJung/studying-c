#include <stdio.h>
#include <stdlib.h>

#define BUFFER_CAPACITY 32

typedef struct {
    int values[BUFFER_CAPACITY];
    size_t start;
    size_t end;
} Circular_Buffer;

void init_buffer(Circular_Buffer *const buf) {
    buf->start = 0;
    buf->end = 0;
}

size_t inc_buf_index(size_t index) {
    index++;
    if (index >= BUFFER_CAPACITY) {
        index = 0;
    }
    return index;
}

size_t dec_buf_index(size_t index) {
    index--;
    if (index >= BUFFER_CAPACITY) {
        index = BUFFER_CAPACITY - 1;
    }
    return index;
}

void print_buf(Circular_Buffer *const buf) {
    printf("Start: %ld, end: %ld\n", buf->start, buf->end);
    size_t index = buf->start;
    while (index != buf->end) {
        printf("[%03ld]: %d\n", index, buf->values[index]);
        index = inc_buf_index(index);
    }
    printf("\n");
}

void push(Circular_Buffer *const buf, const int value) {
    if (buf->start == inc_buf_index(buf->end)) {  // buffer overflow
        return;
    }
    buf->values[buf->end] = value;
    buf->end = inc_buf_index(buf->end);
}

int pop(Circular_Buffer *const buf) {
    if (buf->start == buf->end) {  // buffer is empty
        return -1;
    }
    buf->end = dec_buf_index(buf->end);
    return buf->values[buf->end];
}

// Inserts element at the beginning of the buffer
void unshift(Circular_Buffer *const buf, const int value) {
    size_t index = dec_buf_index(buf->start);
    if (index == buf->end) {  // buffer overflow
        return;
    }
    buf->start = index;
    buf->values[index] = value;
}

// Removes and returns first element of the buffer
int shift(Circular_Buffer *const buf) {
    if (buf->start == buf->end) {  // buffer is empty
        return -1;
    }
    int value = buf->values[buf->start];
    buf->start = inc_buf_index(buf->start);
    return value;
}

int main(void) {
    Circular_Buffer buf;
    init_buffer(&buf);

    printf("Buffer before insertion:\n");
    print_buf(&buf);
    for (int i = 0; i < 16; i++) {
        push(&buf, i);
    }
    printf("Buffer after insertion:\n");
    print_buf(&buf);

    printf("Pop a value from buffer: %d\n", pop(&buf));
    printf("Buffer after pop:\n");
    print_buf(&buf);

    printf("Shift the 1st value from buffer: %d\n", shift(&buf));
    printf("Shift the 2nd value from buffer: %d\n", shift(&buf));
    printf("Buffer after shift:\n");
    print_buf(&buf);

    printf("Buffer after insertion at the beginning:\n");
    for (int i = 0; i < 5; i++) {
        unshift(&buf, 105 - i);
    }
    print_buf(&buf);

    return EXIT_SUCCESS;
}
