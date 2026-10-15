#ifndef PRNG_TINY_MERSENNE_TWISTER_32_H
#define PRNG_TINY_MERSENNE_TWISTER_32_H

void prng_tinymt32_seed(
    /* outputs */
    void* state, /* 16 bytes */
    /*  inputs */
    unsigned long seed /* 32 bits */
);

unsigned long prng_tinymt32_rand(
    /* outputs */
    /* 32 pseudorandom bits */
    /*  inputs */
    void* state
);

#endif /* PRNG_TINY_MERSENNE_TWISTER_32_H */