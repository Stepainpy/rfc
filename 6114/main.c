#include <stdio.h>

#include "cipher_clefia.h"

void memhex(const void* data, size_t size) {
    const unsigned char* byte = data;
    if (!size) return;
    printf("%02x", *byte++); --size;
    while (size --> 0)
        printf(" %02x", *byte++);
}

int main(void) {
    char keys[228], plaintext[16], ciphertext[16];

    cipher_clefia_key_schedule(keys,
        "\xff\xee\xdd\xcc\xbb\xaa\x99\x88\x77\x66\x55\x44\x33\x22\x11\x00", 128);
    cipher_clefia_encrypt(ciphertext, "\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f", keys);
    cipher_clefia_decrypt( plaintext, ciphertext, keys);

    printf("Plaintext : "); memhex( plaintext, 16); putchar('\n');
    printf("Ciphertext: "); memhex(ciphertext, 16); putchar('\n');

    putchar('\n');

    cipher_clefia_key_schedule(keys,
        "\xff\xee\xdd\xcc\xbb\xaa\x99\x88\x77\x66\x55\x44\x33\x22\x11\x00"
        "\xf0\xe0\xd0\xc0\xb0\xa0\x90\x80", 192);
    cipher_clefia_encrypt(ciphertext, "\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f", keys);
    cipher_clefia_decrypt( plaintext, ciphertext, keys);

    printf("Plaintext : "); memhex( plaintext, 16); putchar('\n');
    printf("Ciphertext: "); memhex(ciphertext, 16); putchar('\n');

    putchar('\n');

    cipher_clefia_key_schedule(keys,
        "\xff\xee\xdd\xcc\xbb\xaa\x99\x88\x77\x66\x55\x44\x33\x22\x11\x00"
        "\xf0\xe0\xd0\xc0\xb0\xa0\x90\x80\x70\x60\x50\x40\x30\x20\x10\x00", 256);
    cipher_clefia_encrypt(ciphertext, "\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0a\x0b\x0c\x0d\x0e\x0f", keys);
    cipher_clefia_decrypt( plaintext, ciphertext, keys);

    printf("Plaintext : "); memhex( plaintext, 16); putchar('\n');
    printf("Ciphertext: "); memhex(ciphertext, 16); putchar('\n');

    return 0;
}