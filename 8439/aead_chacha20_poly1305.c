#include "aead_chacha20_poly1305.h"

#include <stdint.h>
#include <string.h>

/* * * * * * * * * * * * * * * * * * * * *
 * Definintion of ChaCha20 stream cipher *
 * * * * * * * * * * * * * * * * * * * * */

typedef uint32_t word_t;

static const word_t chacha20_consts[4] = {
    0x61707865, 0x3320646e, 0x79622d32, 0x6b206574
};

word_t rotl(word_t n, int s) {
    s &= 31; return n << s | n >> (-s & 31);
}

#define chacha20_qround(w, a, b, c, d) do { \
    w[a] += w[b]; w[d] ^= w[a]; w[d] = rotl(w[d], 16); \
    w[c] += w[d]; w[b] ^= w[c]; w[b] = rotl(w[b], 12); \
    w[a] += w[b]; w[d] ^= w[a]; w[d] = rotl(w[d],  8); \
    w[c] += w[d]; w[b] ^= w[c]; w[b] = rotl(w[b],  7); \
} while (0)

void chacha20_block(word_t out[16], word_t counter, const void* key, const void* nonce) {
    word_t init[16], proc[16]; size_t i;

    memcpy(init +  0, chacha20_consts, 16);
    memcpy(init +  4, key, 32);
    memcpy(init + 13, nonce, 12);
    init[12] = counter;

    memcpy(proc, init, 64);

    for (i = 0; i < 10; i++) {
        chacha20_qround(proc, 0, 4,  8, 12);
        chacha20_qround(proc, 1, 5,  9, 13);
        chacha20_qround(proc, 2, 6, 10, 14);
        chacha20_qround(proc, 3, 7, 11, 15);
        chacha20_qround(proc, 0, 5, 10, 15);
        chacha20_qround(proc, 1, 6, 11, 12);
        chacha20_qround(proc, 2, 7,  8, 13);
        chacha20_qround(proc, 3, 4,  9, 14);
    }

    for (i = 0; i < 16; i++)
        out[i] = proc[i] + init[i];
}

/* * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of Poly1305 one-time authenticator *
 * * * * * * * * * * * * * * * * * * * * * * * * */

typedef  int32_t s32_t;
typedef uint32_t u32_t;
typedef  int64_t s64_t;
typedef uint64_t u64_t;

typedef union {
    struct { u32_t lo; s32_t up; } s32;
    struct { u32_t lo, up; } u32;
    s64_t s64;
    u64_t u64;
} register_t;

static const u32_t poly1305_prime[5] = {
    0xfffffffb, 0xffffffff, 0xffffffff, 0xffffffff, 0x3
};

#define poly1305_carw(y, x) do { r.u64 = (u64_t)(x) + r.u32.up; (y) = r.u32.lo; } while (0)

#define poly1305_shlw(y, x, n) do { r.u64 = (u64_t)(x) << (n) | r.u32.up; (y) = r.u32.lo; } while (0)

#define poly1305_addw(z, x, y) do { r.u64 = (u64_t)(x) + (y)           ; (z) = r.u32.lo; } while (0)
#define poly1305_adcw(z, x, y) do { r.u64 = (u64_t)(x) + (y) + r.u32.up; (z) = r.u32.lo; } while (0)

#define poly1305_subw(z, x, y) do { r.s64 = (s64_t)(x) - (y)           ; (z) = r.s32.lo; } while (0)
#define poly1305_sbcw(z, x, y) do { r.s64 = (s64_t)(x) - (y) + r.s32.up; (z) = r.s32.lo; } while (0)

#define poly1305_msfw(y, n, x) do { r.u64 = (u64_t)(n) * (x)           ; (y) = r.u32.lo; } while (0)
#define poly1305_mscw(y, n, x) do { r.u64 = (u64_t)(n) * (x) + r.u32.up; (y) = r.u32.lo; } while (0)

#define poly1305_mafw(z, n, x, y) do { r.u64 = (u64_t)(n) * (x) + (y)           ; (z) = r.u32.lo; } while (0)
#define poly1305_macw(z, n, x, y) do { r.u64 = (u64_t)(n) * (x) + (y) + r.u32.up; (z) = r.u32.lo; } while (0)

u32_t poly1305_shl(u32_t y[5], const u32_t x[5], u32_t n) {
    register_t r = {0};
    poly1305_shlw(y[0], x[0], n);
    poly1305_shlw(y[1], x[1], n);
    poly1305_shlw(y[2], x[2], n);
    poly1305_shlw(y[3], x[3], n);
    poly1305_shlw(y[4], x[4], n);
    return r.u32.up;
}

u32_t poly1305_add(u32_t z[5], const u32_t x[5], const u32_t y[5]) {
    register_t r;
    poly1305_addw(z[0], x[0], y[0]);
    poly1305_adcw(z[1], x[1], y[1]);
    poly1305_adcw(z[2], x[2], y[2]);
    poly1305_adcw(z[3], x[3], y[3]);
    poly1305_adcw(z[4], x[4], y[4]);
    return r.u32.up;
}

s32_t poly1305_sub(u32_t z[5], const u32_t x[5], const u32_t y[5]) {
    register_t r;
    poly1305_subw(z[0], x[0], y[0]);
    poly1305_sbcw(z[1], x[1], y[1]);
    poly1305_sbcw(z[2], x[2], y[2]);
    poly1305_sbcw(z[3], x[3], y[3]);
    poly1305_sbcw(z[4], x[4], y[4]);
    return r.s32.up;
}

void poly1305_mod_prime(u32_t y[5], const u32_t x[10]) {
    u32_t t[5]; size_t i;
    memset(y, 0, 20);
    for (i = 260; i --> 0;) {
        poly1305_shl(y, y, 1);
        y[0] |= x[i / 32] >> (i % 32) & 1;
        if (poly1305_sub(t, poly1305_prime, y))
            poly1305_sub(y, y, poly1305_prime);
    }
}

void poly1305_mul_set(u32_t y[6], u32_t n, const u32_t x[5]) {
    register_t r;
    poly1305_msfw(y[0], n, x[0]);
    poly1305_mscw(y[1], n, x[1]);
    poly1305_mscw(y[2], n, x[2]);
    poly1305_mscw(y[3], n, x[3]);
    poly1305_mscw(y[4], n, x[4]);
    y[5] = r.u32.up;
}

void poly1305_mul_add(u32_t y[6], u32_t n, const u32_t x[5]) {
    register_t r;
    poly1305_mafw(y[0], n, x[0], y[0]);
    poly1305_macw(y[1], n, x[1], y[1]);
    poly1305_macw(y[2], n, x[2], y[2]);
    poly1305_macw(y[3], n, x[3], y[3]);
    poly1305_macw(y[4], n, x[4], y[4]);
    y[5] = r.u32.up;
}

void poly1305_mul(u32_t z[5], const u32_t x[5], const u32_t y[5]) {
    u32_t t[10];

    poly1305_mul_set(t + 0, x[0], y);
    poly1305_mul_add(t + 1, x[1], y);
    poly1305_mul_add(t + 2, x[2], y);
    poly1305_mul_add(t + 3, x[3], y);
    poly1305_mul_add(t + 4, x[4], y);

    poly1305_mod_prime(z, t);
}

/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 * Definition of authenticated encryption with additional data *
 * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

void aead_update_acc(u32_t acc[5], const void* data, const u32_t r[5]) {
    u32_t b[5] = {0};
    memcpy(b, data, 16); b[4] = 1;

    poly1305_add(acc, acc, b);
    poly1305_mul(acc, acc, r);
}

void aead_chacha20_poly1305_encrypt(
    /* outputs */
    void* ciphertext, /* `plainsize` bytes */
    void* tag,        /* 16 bytes */
    /*  inputs */
    const void* plaintext, size_t plainsize,
    const void* aad,       size_t aadsize,
    const void* key,  /* 32 bytes */
    const void* nonce /* 12 bytes */
) {
    u32_t r[5] = {0}, s[5] = {0}, acc[5] = {0};
    word_t gamma[16]; size_t i, c;
    union {
        unsigned char bytes[16];
        word_t words[4];
    } block;

    chacha20_block(gamma, 0, key, nonce);
    memcpy(r, gamma + 0, 16);
    memcpy(s, gamma + 4, 16);
    r[3] &= 0x0ffffffc;
    r[2] &= 0x0ffffffc;
    r[1] &= 0x0ffffffc;
    r[0] &= 0x0fffffff;

    for (i = 0; i < (aadsize + 15) / 16; i++) {
        size_t len = aadsize - 16 * i;
        if (len > 16) len = 16;
        memset(block.bytes, 0, 16);
        memcpy(block.bytes, aad, len);

        aead_update_acc(acc, block.bytes, r);
        aad = (const char*)aad + len;
    }

    for (i = c = 0; i < (plainsize + 15) / 16; i++) {
        size_t len = plainsize - 16 * i;
        if (len > 16) len = 16;
        memcpy(block.bytes, plaintext, len);

        if (i % 4 == 0)
            chacha20_block(gamma, ++c, key, nonce);

        block.words[0] ^= gamma[0 + i % 4 * 4];
        block.words[1] ^= gamma[1 + i % 4 * 4];
        block.words[2] ^= gamma[2 + i % 4 * 4];
        block.words[3] ^= gamma[3 + i % 4 * 4];
        memset(block.bytes + len, 0, 16 - len);
        memcpy(ciphertext, block.bytes, len);

        aead_update_acc(acc, block.bytes, r);
         plaintext = (const char*) plaintext + len;
        ciphertext = (      char*)ciphertext + len;
    }

    memcpy(block.bytes + 0, &  aadsize, 8);
    memcpy(block.bytes + 8, &plainsize, 8);
    aead_update_acc(acc, block.bytes, r);

    poly1305_add(acc, acc, s);

    memcpy(tag, acc, 16);
}

int aead_chacha20_poly1305_decrypt(
    /* outputs */
    void* plaintext, /* `ciphersize` bytes */
    /*  inputs */
    const void* ciphertext, size_t ciphersize,
    const void* aad,        size_t aadsize,
    const void* key,   /* 32 bytes */
    const void* nonce, /* 12 bytes */
    const void* tag    /* 16 bytes */
) {
    u32_t r[5] = {0}, s[5] = {0}, acc[5] = {0};
    word_t gamma[16]; size_t i, c;
    union {
        unsigned char bytes[16];
        word_t words[4];
    } block;

    chacha20_block(gamma, 0, key, nonce);
    memcpy(r, gamma + 0, 16);
    memcpy(s, gamma + 4, 16);
    r[3] &= 0x0ffffffc;
    r[2] &= 0x0ffffffc;
    r[1] &= 0x0ffffffc;
    r[0] &= 0x0fffffff;

    for (i = 0; i < (aadsize + 15) / 16; i++) {
        size_t len = aadsize - 16 * i;
        if (len > 16) len = 16;
        memset(block.bytes, 0, 16);
        memcpy(block.bytes, aad, len);

        aead_update_acc(acc, block.bytes, r);
        aad = (const char*)aad + len;
    }

    for (i = c = 0; i < (ciphersize + 15) / 16; i++) {
        size_t len = ciphersize - 16 * i;
        if (len > 16) len = 16;
        memset(block.bytes, 0, 16);
        memcpy(block.bytes, ciphertext, len);

        if (i % 4 == 0)
            chacha20_block(gamma, ++c, key, nonce);

        aead_update_acc(acc, block.bytes, r);

        block.words[0] ^= gamma[0 + i % 4 * 4];
        block.words[1] ^= gamma[1 + i % 4 * 4];
        block.words[2] ^= gamma[2 + i % 4 * 4];
        block.words[3] ^= gamma[3 + i % 4 * 4];
        memcpy(plaintext, block.bytes, len);

         plaintext = (      char*) plaintext + len;
        ciphertext = (const char*)ciphertext + len;
    }

    memcpy(block.bytes + 0, &   aadsize, 8);
    memcpy(block.bytes + 8, &ciphersize, 8);
    aead_update_acc(acc, block.bytes, r);

    poly1305_add(acc, acc, s);

    return memcmp(acc, tag, 16) == 0;
}