#include "prng_tinymt32.h"

#include <stdint.h>

typedef struct { uint32_t s[4]; } tinymt32_t;

void tinymt32_step(tinymt32_t* tmt) {
    uint32_t x, y;

    y = tmt->s[3];
    x = (tmt->s[0] & UINT32_C(0x7fffffff)) ^ tmt->s[1] ^ tmt->s[2];
    x ^= x << 1;
    y ^= (y >> 1) ^ x;

    tmt->s[0] = tmt->s[1];
    tmt->s[1] = tmt->s[2];
    tmt->s[2] = x ^ (y << 10);
    tmt->s[3] = y;

    if (y & 1) {
        tmt->s[1] ^= UINT32_C(0x8f7011ee);
        tmt->s[2] ^= UINT32_C(0xfc78ff1f);
    }
}

void prng_tinymt32_seed(void* state, unsigned long seed) {
    tinymt32_t* tmt = state; int i;

    tmt->s[0] = seed;
    tmt->s[1] = UINT32_C(0x8f7011ee);
    tmt->s[2] = UINT32_C(0xfc78ff1f);
    tmt->s[3] = UINT32_C(0x3793fdff);

    for (i = 1; i < 8; i++)
        tmt->s[i & 3] ^= i + UINT32_C(1812433253) *
            (tmt->s[(i - 1) & 3] ^ tmt->s[(i - 1) & 3] >> 30);

    for (i = 0; i < 8; i++)
        tinymt32_step(tmt);
}

unsigned long prng_tinymt32_rand(void* state) {
    tinymt32_t* tmt = state;
    uint32_t t0, t1;

    tinymt32_step(tmt);

    t0 = tmt->s[3];
    t1 = tmt->s[0] + (tmt->s[2] >> 8);
    t0 ^= t1;
    if (t1 & 1)
        t0 ^= UINT32_C(0x3793fdff);

    return t0;
}