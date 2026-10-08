#ifndef AEAD_CHACHA20_POLY1305_H
#define AEAD_CHACHA20_POLY1305_H

#include <stddef.h>

void aead_chacha20_poly1305_encrypt(
    /* outputs */
    void* ciphertext, /* `plainsize` bytes */
    void* tag,        /* 16 bytes */
    /*  inputs */
    const void* plaintext, size_t plainsize,
    const void* aad,       size_t aadsize,
    const void* key,  /* 32 bytes */
    const void* nonce /* 12 bytes */
);

int aead_chacha20_poly1305_decrypt(
    /* outputs */
    void* plaintext, /* `ciphersize` bytes */
    /*  inputs */
    const void* ciphertext, size_t ciphersize,
    const void* aad,        size_t aadsize,
    const void* key,   /* 32 bytes */
    const void* nonce, /* 12 bytes */
    const void* tag    /* 16 bytes */
);

#endif /* AEAD_CHACHA20_POLY1305_H */