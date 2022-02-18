#include "../../lib/integration.h"

#include <stdio.h>
#include <stdlib.h>

double f(const double x) {
    return (x * x);
}

int main(void) {
    const double start = 0.0;
    const double end = 1.0;
    const size_t num_of_segments = 10;
    printf("Integration from %g to %g results in: %g (Trapezoidal rule)\n", start, end, trapezoidal_rule(f, start, end, num_of_segments));
    printf("Integration from %g to %g results in: %g (Simpson rule)\n", start, end, simpson_rule(f, start, end, num_of_segments));
    printf("Integration from %g to %g results in: %g (Simpson rule 3/8)\n", start, end, simpson_rule_38(f, start, end, num_of_segments));
    return EXIT_SUCCESS;
}
