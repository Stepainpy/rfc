#include <stdio.h>

#include "prng_tinymt32.h"

int main(void) {
    char s[16]; int i;

    prng_tinymt32_seed(s, 1);
    for (i = 1; i <= 50; i++) {
        printf("%10lu ", prng_tinymt32_rand(s));
        if (i % 5 == 0) putchar('\n');
    }

    return 0;
}