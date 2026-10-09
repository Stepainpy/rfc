#ifndef TINY_MERSENNE_TWISTER_32_H
#define TINY_MERSENNE_TWISTER_32_H

#include <stdint.h>

typedef struct { uint32_t s[4]; } tinymt32_t;

void tinymt32_seed(tinymt32_t* state, uint32_t seed);

uint32_t tinymt32_rand(tinymt32_t* state);

#endif /* TINY_MERSENNE_TWISTER_32_H */