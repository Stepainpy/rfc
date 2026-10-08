#include <stdio.h>
#include <ctype.h>
#include <string.h>

#include "aead_chacha20_poly1305.h"

void chacha20_block(void* o, unsigned c, const void* k, const void* n);

void poly1305_add(void* z, const void* x, const void* y);
void poly1305_mul(void* z, const void* x, const void* y);

void memhex(const void* data, size_t size) {
    const unsigned char* byte = data;
    if (!size) return;
    printf("%02x", *byte++); --size;
    while (size --> 0)
        printf(" %02x", *byte++);
}

void hexdump(const void* ptr, size_t size) {
    const unsigned char* data = ptr;
    size_t offset = 0, i, j;
    for (i = 0; i < (size + 15) / 16; i++) {
        size_t len = size - offset < 16 ? size - offset : 16;
        printf("%03u ", (unsigned)offset);
        memhex(data + offset, len);
        for (j = len; j < 16; j++) printf("   ");
        printf("  ");
        for (j = 0; j < len; j++)
            printf("%c", isprint(data[offset + j]) ? data[offset + j] : '.');
        putchar('\n');
        offset += len;
    }
}

int main_chacha(void) {
    unsigned gamma[16];

    chacha20_block(gamma, 1,
        "\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f"
        "\x10\x11\x12\x13\x14\x15\x16\x17\x18\x19\x1a\x1b\x1c\x1d\x1e\x1f",
        "\x00\x00\x00\x09\x00\x00\x00\x4a\x00\x00\x00\x00"
    );
    hexdump(gamma, 64);

    chacha20_block(gamma, 0,
        "\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f"
        "\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f",
        "\x00\x00\x00\x00\x00\x01\x02\x03\x04\x05\x06\x07"
    );
    hexdump(gamma, 64);

    return 0;
}

int main_poly_base(void) {
    FILE* tv = fopen("tv.txt", "r"); int i, e;
    if (!tv) return 1;

    for (i = e = 0; i < 100; i++) {
        unsigned a[5], b[5], r[5], s[5], p[5];
        unsigned S[5], P[5];

        fscanf(tv, "%1x%8x%8x%8x%8x", a+4, a+3, a+2, a+1, a); fgetc(tv);
        fscanf(tv, "%1x%8x%8x%8x%8x", b+4, b+3, b+2, b+1, b); fgetc(tv);
        fscanf(tv, "%1x%8x%8x%8x%8x", r+4, r+3, r+2, r+1, r); fgetc(tv);
        fscanf(tv, "%1x%8x%8x%8x%8x", s+4, s+3, s+2, s+1, s); fgetc(tv);
        fscanf(tv, "%1x%8x%8x%8x%8x", p+4, p+3, p+2, p+1, p); fgetc(tv);

        poly1305_add(S, a, b);
        poly1305_mul(P, S, r);

        if (memcmp(S, s, 20) != 0) {
            printf("%2i) wrong sum\n", i); ++e;
            printf("X: "); memhex(S, 20); putchar('\n');
            printf("O: "); memhex(s, 20); putchar('\n');
        }
        if (memcmp(P, p, 20) != 0) {
            printf("%2i) wrong product\n", i); ++e;
            printf("X: "); memhex(P, 20); putchar('\n');
            printf("O: "); memhex(p, 20); putchar('\n');
        }
        if (e >= 10) break;
     }

    fclose(tv);
    return 0;
}

int main_poly_loop(void) {
    unsigned acc[5] = {0}, msg[15] = {0}, i;
    const unsigned r[5] = { 0x08bed685, 0x036d5554, 0x0e52447c, 0x0806d540, 0 };
    const unsigned s[5] = { 0x8a800301, 0xfdb20dfb, 0xaff6bf4a, 0x1bf54941, 0 };

    memcpy(msg +  0, "Cryptographic Fo" "\1", 17);
    memcpy(msg +  5, "rum Research Gro" "\1", 17);
    memcpy(msg + 10, "up"               "\1",  3);

    for (i = 0; i < 3; i++) {
        poly1305_add(acc, acc, msg + 5*i);
        poly1305_mul(acc, acc, r);
    }
    poly1305_add(acc, acc, s);

    hexdump(acc, 16);
    return 0;
}

int main(void) {
    const char* plaintext =
        "Ladies and Gentlemen of the class of '99: If I could offer "
        "you only one tip for the future, sunscreen would be it.";
    const char* aad =
        "\x50\x51\x52\x53\xc0\xc1\xc2\xc3\xc4\xc5\xc6\xc7";

    const char* key =
        "\x80\x81\x82\x83\x84\x85\x86\x87\x88\x89\x8a\x8b\x8c\x8d\x8e\x8f"
        "\x90\x91\x92\x93\x94\x95\x96\x97\x98\x99\x9a\x9b\x9c\x9d\x9e\x9f";
    const char* nonce = "\x07\x00\x00\x00\x40\x41\x42\x43\x44\x45\x46\x47";

    char ciphertext[114], tag[16];
    char  plaindata[114]; int verify;

    aead_chacha20_poly1305_encrypt(ciphertext, tag, plaintext, 114, aad, 12, key, nonce);

    puts("Ciphertext:");
    hexdump(ciphertext, 114);

    puts("\nTag:");
    hexdump(tag, 16);

    verify = aead_chacha20_poly1305_decrypt(plaindata, ciphertext, 114, aad, 12, key, nonce, tag);

    puts("\nPlaintext:");
    hexdump(plaindata, 114);

    printf("\nTag: %s\n", verify ? "valid" : "invalid");

    return 0;
}