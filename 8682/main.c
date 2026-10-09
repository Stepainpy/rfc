#include <stdio.h>

#include "tinymt32.h"

int main(void) {
    tinymt32_t s; int i;

    tinymt32_seed(&s, 1);
    for (i = 1; i <= 50; i++) {
        printf("%10u ", tinymt32_rand(&s));
        if (i % 5 == 0) putchar('\n');
    }

    return 0;
}