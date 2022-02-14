#include <stdio.h>
#include <stdlib.h>

typedef struct Node Node;
struct Node {
    int value;
    Node *previous;
    Node *next;
};

void print_list(const Node *const head) {
    size_t index = 0;
    const Node *current = head;
    printf("Index - Address - Value\n");
    while (current != NULL) {
        printf("%03lu - %p - %d\n", index, (void *)current, current->value);
        current = current->next;
        index++;
    }
    printf("\n");
}

Node *insert_before(Node *const next, const int value) {
    Node *node = malloc(sizeof(Node));
    if (node != NULL) {
        node->value = value;
        node->next = next;
        if (next != NULL) {
            Node *previous = next->previous;
            node->previous = previous;
            next->previous = node;
            if (previous != NULL) {
                previous->next = node;
            }
        } else {
            node->previous = NULL;
        }
    }
    return node;
}

Node *insert_after(Node *const previous, const int value) {
    Node *node = malloc(sizeof(Node));
    if (node != NULL) {
        node->value = value;
        node->previous = previous;
        if (previous != NULL) {
            Node *next = previous->next;
            node->next = next;
            previous->next = node;
            if (next != NULL) {
                next->previous = node;
            }
        } else {
            node->next = NULL;
        }
    }
    return node;
}

Node *find_head(Node *const node) {
    // Sequential search
    Node *current = node;
    while (current != NULL && current->previous != NULL) {
        current = current->previous;
    }
    return current;
}

Node *find_tail(Node *const node) {
    // Sequential search
    Node *current = node;
    while (current != NULL && current->next != NULL) {
        current = current->next;
    }
    return current;
}

Node *find_in_list(Node *const head, const int value) {
    // Sequential search
    Node *current = head;
    while (current != NULL) {
        if (current->value == value) {
            return current;
        }
        current = current->next;
    }
    return NULL;
}

Node *reverse_find_in_list(Node *const tail, const int value) {
    // Sequential search
    Node *current = tail;
    while (current != NULL) {
        if (current->value == value) {
            return current;
        }
        current = current->previous;
    }
    return NULL;
}

Node *insert_at_beginning(Node *const node, const int value) {
    Node *head = find_head(node);
    return insert_before(head, value);
}

Node *insert_at_end(Node *const node, const int value) {
    Node *tail = find_tail(node);
    return insert_after(tail, value);
}

void delete_in_list(Node *const node) {
    if (node == NULL) {
        return;
    }
    Node *previous = node->previous;
    Node *next = node->next;
    if (previous != NULL) {
        previous->next = next;
    }
    if (next != NULL) {
        next->previous = previous;
    }
    free(node);
}

void free_list(Node **const head) {
    if (head == NULL || *head == NULL) {
        return;
    }
    Node *current = *head;
    while (current != NULL) {
        Node *next = current->next;
        free(current);
        current = next;
    }
    *head = NULL;
}

int main(void) {
    Node *second = NULL;
    printf("Double linked list before insertion:\n");
    print_list(second);
    second = insert_before(NULL, 2);
    Node *first = insert_before(second, 1);
    Node *third = insert_after(second, 3);
    Node *fourth = insert_after(third, 4);
    insert_after(fourth, 5);
    printf("Double linked list after insertion:\n");
    print_list(first);

    printf("find 8 in list results in address %p\n", (void *)find_in_list(second, 8));
    printf("find NULL in list results in address %p\n", (void *)find_in_list(NULL, 3));
    Node *found = find_in_list(second, 3);
    printf("find 3 in list results in address %p with value %d\n\n", (void *)found, found->value);

    printf("find head from NULL pointer results in address %p\n", (void *)find_head(NULL));
    Node *head = find_head(fourth);
    printf("find head in list results in address %p with value %d\n\n", (void *)head, head->value);

    printf("find tail from NULL pointer results in address %p\n", (void *)find_tail(NULL));
    Node *tail = find_tail(third);
    printf("find tail in list results in address %p with value %d\n\n", (void *)tail, tail->value);

    delete_in_list(NULL);
    delete_in_list(found);
    printf("\nDouble linked list after deletion of element with value 3:\n");
    print_list(head);

    printf("reverse find 8 in list results in address %p\n", (void *)reverse_find_in_list(tail, 8));
    printf("reverse find NULL in list results in address %p\n", (void *)reverse_find_in_list(NULL, 3));
    found = reverse_find_in_list(tail, 4);
    printf("reverse find 4 in list results in address %p with value %d\n\n", (void *)found, found->value);

    free_list(&head->next->next);
    printf("Double linked list after deallocation of two last nodes:\n");
    print_list(head);
    free_list(NULL);

    found = reverse_find_in_list(find_tail(first), 1);
    printf("reverse find 1 in list results in address %p with value %d\n", (void *)found, found->value);
    found = find_in_list(find_head(second), 2);
    printf("find 2 in list results in address %p with value %d\n\n", (void *)found, found->value);

    printf("Double linked list after deallocation of full list:\n");
    free_list(&head);
    print_list(head);

    printf("Double linked list after insertion of new elements:\n");
    second = insert_before(NULL, 2);
    insert_at_beginning(second, 1);
    insert_at_end(second, 3);
    insert_at_end(second, 4);
    head = find_head(second);
    print_list(head);
    free_list(&head);
    return EXIT_SUCCESS;
}
