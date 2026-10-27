#ifndef EDDSA_ED25519_H
#define EDDSA_ED25519_H

#include <stddef.h>

void eddsa_ed25519_key_generation(
    /* outputs */
    void* private_key, /* 32 bytes */
    void* public_key,  /* 32 bytes */
    /*  inputs */
    const void* random /* 32 bytes */
);

int eddsa_ed25519_sign(
    /* outputs */
    void* signature, /* 64 bytes */
    /*  inputs */
    const void* message, size_t msgsize,
    const void* context, size_t ctxsize, /* nullable or 0-255 bytes */
    const void* private_key, /* 32 bytes */
    const void* public_key   /* 32 bytes */
);

int eddsa_ed25519_verify(
    /* outputs */
    /*  inputs */
    const void* signature, /* 64 bytes */
    const void* message, size_t msgsize,
    const void* context, size_t ctxsize, /* nullable or 0-255 bytes */
    const void* public_key /* 32 bytes */
);

#endif /* EDDSA_ED25519_H */