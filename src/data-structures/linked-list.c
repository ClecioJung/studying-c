#include <stdio.h>
#include <stdlib.h>

typedef struct Node Node;
struct Node {
    int value;
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

Node *insert_new_head(Node **const head, const int value) {
    Node *node = malloc(sizeof(Node));
    if (node != NULL) {
        node->value = value;
        node->next = *head;
        *head = node;
    }
    return node;
}

Node *insert_after(Node **const previous, const int value) {
    Node *node = malloc(sizeof(Node));
    if (node != NULL) {
        node->value = value;
        if (*previous != NULL) {
            node->next = (*previous)->next;
            (*previous)->next = node;
        } else {
            node->next = NULL;
            *previous = node;
        }
    }
    return node;
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

void delete_next(Node *const previous) {
    if (previous == NULL) {
        return;
    }
    Node *node = previous->next;
    if (node == NULL) {
        return;
    }
    previous->next = node->next;
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
    Node *head = NULL;
    printf("Linked list before insertion:\n");
    print_list(head);
    Node *last = insert_after(&head, 3);
    insert_new_head(&head, 2);
    insert_new_head(&head, 1);
    last = insert_after(&last, 4);
    insert_after(&last, 5);
    printf("Linked list after insertion:\n");
    print_list(head);

    Node *found = find_in_list(head, 8);
    printf("find 8 in list results in address %p\n", (void *)found);
    found = find_in_list(NULL, 3);
    printf("find NULL in list results in address %p\n", (void *)found);
    found = find_in_list(head, 3);
    printf("find 3 in list results in address %p with value %d\n", (void *)found, found->value);

    delete_next(NULL);
    delete_next(found);
    printf("\nLinked list after deletion of element with value 4:\n");
    print_list(head);

    free_list(&head->next->next);
    printf("Linked list after deallocation of two last nodes:\n");
    print_list(head);

    free_list(NULL);
    free_list(&head);
    printf("Linked list after deallocation of full list:\n");
    print_list(head);

    printf("Linked list after insertion of new head:\n");
    insert_new_head(&head, 1);
    print_list(head);
    free_list(&head);
    return EXIT_SUCCESS;
}
