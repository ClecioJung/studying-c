#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct Node Node;
struct Node {
    char value;
    Node *left;
    Node *right;
};

void print_node(const Node *const node, const unsigned int level) {
    if (node == NULL) {
        return;
    }
    printf("%c\n", node->value);
    if (node->left) {
        printf("%*sLEFT:  ", level * 2, "");
        print_node(node->left, level + 1);
    }
    if (node->right) {
        printf("%*sRIGHT: ", level * 2, "");
        print_node(node->right, level + 1);
    }
}

void print_tree(const Node *const root) {
    printf("HEAD:  ");
    print_node(root, 0);
}

void free_tree(Node **const root) {
    if (root == NULL || *root == NULL) {
        return;
    }
    free_tree(&((*root)->left));
    free_tree(&((*root)->right));
    free(*root);
    *root = NULL;
}

Node *insert_new_root(Node **const root, const char value) {
    Node *node = malloc(sizeof(Node));
    if (node != NULL) {
        node->value = value;
        node->left = *root;
        node->right = NULL;
        *root = node;
    }
    return node;
}

Node *insert_leaf_right(Node **const root, const char value) {
    Node *node = malloc(sizeof(Node));
    if (node == NULL) {
        return NULL;
    }
    node->value = value;
    node->right = NULL;
    node->left = NULL;
    if (*root == NULL) {
        *root = node;
    } else {
        Node *previous = *root;
        while (previous->right != NULL) {
            previous = previous->right;
        }
        previous->right = node;
    }
    return node;
}

Node *insert_node_after(Node *const previous, const char value) {
    if (previous == NULL) {
        return NULL;
    }
    Node *node = malloc(sizeof(Node));
    if (node == NULL) {
        return NULL;
    }
    Node *after = previous->right;
    node->value = value;
    node->right = NULL;
    node->left = after;
    previous->right = node;
    return node;
}

int precedence(const char op) {
    static const char ops[] = {'^', '%', '/', '*', '-', '+'};
    static const int num_ops = sizeof(ops) / sizeof(ops[0]);
    for (int index = 0; index < num_ops; index++) {
        if (op == ops[index]) {
            return index;
        }
    }
    return -1;
}

Node *find_previous_right(Node *const root, const char value) {
    if (root == NULL) {
        return NULL;
    }
    Node *previous = root;
    const char prec = precedence(value);
    while (previous->right != NULL) {
        if (precedence(previous->right->value) <= prec) {
            break;
        }
        previous = previous->right;
    }
    return previous;
}

Node *insert_new_op(Node **const root, const char value) {
    if (root == NULL) {
        return NULL;
    }
    Node *previous = find_previous_right(*root, value);
    if (previous != NULL && precedence(previous->value) > precedence(value)) {
        return insert_node_after(previous, value);
    } else {
        return insert_new_root(root, value);
    }
}

Node *build_ast(const char *const expr) {
    Node *root = NULL;
    const int len = strlen(expr);
    for (int i = 0; i < len; i++) {
        const char c = expr[i];
        if (isalnum(c)) {
            insert_leaf_right(&root, c);
        } else {
            insert_new_op(&root, c);
        }
    }
    return root;
}

int power(int base, int expoent) {
    int result = 1;
    while (expoent != 0) {
        result *= base;
        expoent--;
    }
    return result;
}

int eval(Node *const node) {
    if (node == NULL) {
        return 0;
    }
    const char c = node->value;
    switch (c) {
        case '+':
            return eval(node->left) + eval(node->right);
        case '-':
            return eval(node->left) - eval(node->right);
        case '*':
            return eval(node->left) * eval(node->right);
        case '/':
            return eval(node->left) / eval(node->right);
        case '%':
            return eval(node->left) % eval(node->right);
        case '^':
            return power(eval(node->left), eval(node->right));
        default:
            return (c - '0');
    }
}

void solve_expr(const char *const expr) {
    Node *root = build_ast(expr);
    printf("Expr %s results in: %d\n", expr, eval(root));
    printf("Binary tree:\n");
    print_tree(root);
    printf("\n");
    free_tree(&root);
}

/* Simple expression parser developed as an example to study binary trees
   Limitations:
     - numbers are single character
     - it performs the operation with integer numbers
     - accept only the following operations: +, -, *, /, % (modulus) and ^ (power)
*/
int main(void) {
    solve_expr("5*2+3");
    solve_expr("3+5*2");
    solve_expr("5*3+1*9+9/3");
    solve_expr("6/2/3-7*2/2");
    solve_expr("5*3-2*2+2^1");
    solve_expr("5*3+1*9+9/3-6/2/3-7*2/2");
    solve_expr("9%2+2-1");
    solve_expr("2");
    return EXIT_SUCCESS;
}
