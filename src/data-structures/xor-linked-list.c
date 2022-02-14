#include <stdio.h>
#include <stdlib.h>

typedef struct Node Node;
struct Node {
    int value;
    Node *xor_addr;
};

Node *next_node(const Node *const previous, const Node *const node) {
    if (node == NULL) {
        return NULL;
    }
    return (Node *)((size_t)previous ^ (size_t)node->xor_addr);
}

Node *previous_node(const Node *const node, const Node *const next) {
    if (node == NULL) {
        return NULL;
    }
    return (Node *)((size_t)next ^ (size_t)node->xor_addr);
}

void set_address(Node **const node, Node *const previous, Node *const next) {
    (*node)->xor_addr = (Node *)((size_t)previous ^ (size_t)next);
    if (next != NULL) {
        Node *after_next = next_node(previous, next);
        next->xor_addr = (Node *)((size_t)(*node) ^ (size_t)after_next);
    }
    if (previous != NULL) {
        Node *before_previous = previous_node(previous, next);
        previous->xor_addr = (Node *)((size_t)(*node) ^ (size_t)before_previous);
    }
}

void print_list(const Node *const head) {
    size_t index = 0;
    const Node *previous = NULL;
    const Node *current = head;
    printf("Index - Address - Value\n");
    while (current != NULL) {
        printf("%03lu - %p - %d\n", index, (void *)current, current->value);
        const Node *next = next_node(previous, current);
        previous = current;
        current = next;
        index++;
    }
    printf("\n");
}

Node *find_in_list(Node *const head, const int value, Node **const previous) {
    // Sequential search
    *previous = NULL;
    Node *current = head;
    while (current != NULL) {
        if (current->value == value) {
            return current;
        }
        Node *next = next_node(*previous, current);
        *previous = current;
        current = next;
    }
    return NULL;
}

Node *reverse_find_in_list(Node *const tail, const int value, Node **const next) {
    // Sequential search
    *next = NULL;
    Node *current = tail;
    while (current != NULL) {
        if (current->value == value) {
            return current;
        }
        Node *previous = previous_node(current, *next);
        *next = current;
        current = previous;
    }
    return NULL;
}

Node *find_head(Node *const tail) {
    // Sequential search
    Node *next = tail;
    Node *current = previous_node(tail, NULL);
    while (current != NULL) {
        Node *previous = previous_node(current, next);
        next = current;
        current = previous;
    }
    return next;
}

Node *find_tail(Node *const head) {
    // Sequential search
    Node *previous = head;
    Node *current = next_node(NULL, head);
    while (current != NULL) {
        Node *next = next_node(previous, current);
        previous = current;
        current = next;
    }
    return previous;
}

Node *insert_at_beginning(Node **const head, const int value) {
    if (head == NULL) {
        return NULL;
    }
    Node *node = malloc(sizeof(Node));
    if (node != NULL) {
        node->value = value;
        set_address(&node, NULL, *head);
        *head = node;
    }
    return node;
}

Node *insert_new_tail(Node *const tail, const int value) {
    if (tail == NULL) {
        return NULL;
    }
    Node *node = malloc(sizeof(Node));
    if (node != NULL) {
        node->value = value;
        set_address(&node, tail, NULL);
    }
    return node;
}

Node *insert_at_end(Node **const head, const int value) {
    if (head == NULL) {
        return NULL;
    }
    Node *tail = find_tail(*head);
    return insert_new_tail(tail, value);
}

Node *insert_between(Node *const previous, Node *const next, const int value) {
    if (previous == NULL && next == NULL) {
        return NULL;
    }
    Node *node = malloc(sizeof(Node));
    if (node != NULL) {
        node->value = value;
        set_address(&node, previous, next);
    }
    return node;
}

void delete_in_list(Node *const node, Node *const previous) {
    if (node == NULL) {
        return;
    }
    Node *next = next_node(previous, node);
    if (next != NULL) {
        Node *after_next = next_node(node, next);
        next->xor_addr = (Node *)((size_t)(previous) ^ (size_t)after_next);
    }
    if (previous != NULL) {
        Node *before_previous = previous_node(previous, node);
        previous->xor_addr = (Node *)((size_t)(next) ^ (size_t)before_previous);
    }
    free(node);
}

void free_list(Node **const head) {
    if (head == NULL || *head == NULL) {
        return;
    }
    const Node *previous = NULL;
    Node *current = *head;
    while (current != NULL) {
        Node *next = next_node(previous, current);
        free(current);
        previous = current;
        current = next;
    }
    *head = NULL;
}

int main(void) {
    Node *head = NULL;
    printf("XOR linked list before insertion:\n");
    print_list(head);
    insert_at_beginning(&head, 2);
    insert_at_beginning(&head, 1);
    insert_at_end(&head, 3);
    insert_at_end(&head, 4);
    insert_at_end(&head, 5);
    printf("XOR linked list after insertion:\n");
    print_list(head);

    Node *previous;
    printf("find 8 in list results in address %p\n", (void *)find_in_list(head, 8, &previous));
    printf("find NULL in list results in address %p\n", (void *)find_in_list(NULL, 3, &previous));
    Node *found = find_in_list(head, 3, &previous);
    printf("find 3 in list results in address %p with value %d\n\n", (void *)found, found->value);

    printf("find tail from NULL pointer results in address %p\n", (void *)find_tail(NULL));
    Node *tail = find_tail(head);
    printf("find tail in list results in address %p with value %d\n\n", (void *)tail, tail->value);

    {
        printf("find head from NULL pointer results in address %p\n", (void *)find_head(NULL));
        Node *_head = find_head(tail);
        printf("find head in list results in address %p with value %d\n\n", (void *)_head, _head->value);
    }

    delete_in_list(NULL, NULL);
    delete_in_list(found, previous);
    printf("Double linked list after deletion of element with value 3:\n");
    print_list(head);

    Node *next;
    printf("reverse find 8 in list results in address %p\n", (void *)reverse_find_in_list(tail, 8, &next));
    printf("reverse find NULL in list results in address %p\n", (void *)reverse_find_in_list(NULL, 3, &next));
    found = reverse_find_in_list(tail, 4, &next);
    printf("reverse find 4 in list results in address %p with value %d\n\n", (void *)found, found->value);

    {
        insert_between(found, next, 9);
        printf("Double linked list after insertion of element with value 9:\n");
        print_list(head);
    }

    {
        insert_new_tail(tail, 10);
        printf("Double linked list after insertion of new tail with value 10:\n");
        print_list(head);
    }

    printf("Double linked list after deallocation of full list:\n");
    free_list(NULL);
    free_list(&head);
    print_list(head);

    return EXIT_SUCCESS;
}
