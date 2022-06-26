#include <stdio.h>
#include <stdlib.h>

#include "../../lib/linear-systems.h"

// Remember to free the returned vector after calling this function!
Vector compute_tridiagonal_system(const Vector t, const Vector r, const Vector d, const Vector x) {
    if ((t.len == 0) || (t.len != r.len) || (t.len != d.len) || (t.len != x.len)) {
        return (Vector){0};
    }
    Vector result = vector_alloc(t.len);
    if (t.len == 1) {
        result.data[0] = r.data[0] * x.data[0];
    } else {
        result.data[0] = r.data[0] * x.data[0] + d.data[0] * x.data[1];
        for (size_t i = 1; (i + 1) < t.len; i++) {
            result.data[i] = t.data[i] * x.data[i - 1] + r.data[i] * x.data[i] + d.data[i] * x.data[i + 1];
        }
        result.data[t.len - 1] = t.data[t.len - 1] * x.data[t.len - 2] + r.data[t.len - 1] * x.data[t.len - 1];
    }
    return result;
}

int main(void) {
    Vector t = vector_alloc(5);
    Vector r = vector_alloc(5);
    Vector d = vector_alloc(5);
    Vector b = vector_alloc(5);
    Vector x = (Vector){0};
    t.data[0] = 0.0;
    t.data[1] = 1.0;
    t.data[2] = 1.0;
    t.data[3] = -1.0;
    t.data[4] = -1.0;
    r.data[0] = 1.0;
    r.data[1] = 1.0;
    r.data[2] = -1.0;
    r.data[3] = 1.0;
    r.data[4] = 2.0;
    d.data[0] = -1.0;
    d.data[1] = -1.0;
    d.data[2] = 1.0;
    d.data[3] = 1.0;
    d.data[4] = 0.0;
    b.data[0] = 0.0;
    b.data[1] = 1.0;
    b.data[2] = 2.0;
    b.data[3] = -1.0;
    b.data[4] = -2.0;
    printf("\nSolving a triduagonal linear system\n");
    printf("Vector t:\n");
    vector_print(t);
    printf("Vector r:\n");
    vector_print(r);
    printf("Vector d:\n");
    vector_print(d);
    printf("Vector b:\n");
    vector_print(b);
    if (tridiagonal_solving(t, r, d, b, &x) != Solvable) {
        fprintf(stderr, "Couldn't solve the tridiagonal system!\n\n");
        return EXIT_FAILURE;
    }
    printf("The solution is:\n");
    vector_print(x);
    Vector result = compute_tridiagonal_system(t, r, d, x);
    if (vector_are_equal(result, b)) {
        printf("This is the correct solution!\n");
    } else {
        fprintf(stderr, "The agorithm didn't found the correct solution!\n");
        return EXIT_FAILURE;
    }
    vector_dealloc(&t);
    vector_dealloc(&r);
    vector_dealloc(&d);
    vector_dealloc(&b);
    vector_dealloc(&x);
    vector_dealloc(&result);
    return EXIT_SUCCESS;
}
