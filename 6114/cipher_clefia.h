#ifndef CIPHER_CLEFIA_H
#define CIPHER_CLEFIA_H

int cipher_clefia_key_schedule(
    /* outputs */
    void* key_schedule, /* 164/196/228 bytes */
    /*  inputs */
    const void* key, /* 16/24/32 bytes */
    int bits         /* 128/192/256 */
);

void cipher_clefia_encrypt(
    /* outputs */
    void* ciphertext, /* 16 bytes */
    /*  inputs */
    const void* plaintext,   /*          16 bytes */
    const void* key_schedule /* 164/196/228 bytes */
);

void cipher_clefia_decrypt(
    /* outputs */
    void* plaintext, /* 16 bytes */
    /*  inputs */
    const void* ciphertext,  /*          16 bytes */
    const void* key_schedule /* 164/196/228 bytes */
);

#endif /* CIPHER_CLEFIA_H */