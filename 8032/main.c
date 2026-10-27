#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include "eddsa_ed25519.h"

uint64_t arith_add(uint64_t z[4], const uint64_t x[4], const uint64_t y[4]);
uint64_t arith_sub(uint64_t z[4], const uint64_t x[4], const uint64_t y[4]);
void     arith_mul(uint64_t z[8], const uint64_t x[4], const uint64_t y[4]);

void prime_add(uint64_t z[4], const uint64_t x[4], const uint64_t y[4]);
void prime_sub(uint64_t z[4], const uint64_t x[4], const uint64_t y[4]);
void prime_mul(uint64_t z[4], const uint64_t x[4], const uint64_t y[4]);

void order_add(uint64_t z[4], const uint64_t x[4], const uint64_t y[4]);
void order_mul(uint64_t z[4], const uint64_t x[4], const uint64_t y[4]);

void prime_exp_by_2_252_3(uint64_t y[4], const uint64_t x[4]);
void prime_inv(uint64_t y[4], const uint64_t x[4]);

void prime_mod(uint64_t x[4]);
void order_mod(uint64_t x[4]);

void point_pack(void* out, const void* p);
void point_unpack(void* p, const void* in);

void shahf_begin(void* ctx);
void shahf_update(void* ctx, const void* data, size_t size);
void shahf_end(void* ctx, void* hash);

int test_add(void) {
    FILE* file = fopen("tests/add.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], z[4], Z[4], c, C;

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);
        c = fgetc(file) - '0';
        fscanf(file, "%16lx%16lx%16lx%16lx", z+3, z+2, z+1, z+0); fgetc(file);

        C = arith_add(Z, x, y);

        if (C != c || memcmp(Z, z, 32) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : %01lu%016lx%016lx%016lx%016lx\n", C, Z[3], Z[2], Z[1], Z[0]);
            fprintf(stderr, "  correct: %01lu%016lx%016lx%016lx%016lx\n", c, z[3], z[2], z[1], z[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_sub(void) {
    FILE* file = fopen("tests/sub.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], z[4], Z[4], c, C;

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);
        c = fgetc(file) - '0';
        fscanf(file, "%16lx%16lx%16lx%16lx", z+3, z+2, z+1, z+0); fgetc(file);

        C = -arith_sub(Z, x, y);

        if (C != c || memcmp(Z, z, 32) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : %01lu%016lx%016lx%016lx%016lx\n", C, Z[3], Z[2], Z[1], Z[0]);
            fprintf(stderr, "  correct: %01lu%016lx%016lx%016lx%016lx\n", c, z[3], z[2], z[1], z[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_mul(void) {
    FILE* file = fopen("tests/mul.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], z[8], Z[8];

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx%16lx%16lx%16lx%16lx",
            z+7, z+6, z+5, z+4, z+3, z+2, z+1, z+0); fgetc(file);

        arith_mul(Z, x, y);

        if (memcmp(Z, z, 64) != 0) {
            size_t k;
            fprintf(stderr, "Case #%2lu   ", i);
            for (k = 128; k --> 0;) { putchar((Z[k/16] >> (k%16*4) & 15) != (z[k/16] >> (k%16*4) & 15) ? '_' : ' '); } putchar('\n');
            fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx%016lx%016lx%016lx%016lx\n", Z[7], Z[6], Z[5], Z[4], Z[3], Z[2], Z[1], Z[0]);
            fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx%016lx%016lx%016lx%016lx\n", z[7], z[6], z[5], z[4], z[3], z[2], z[1], z[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_addp(void) {
    FILE* file = fopen("tests/addp.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], z[4], Z[4];

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", z+3, z+2, z+1, z+0); fgetc(file);

        prime_add(Z, x, y);

        if (memcmp(Z, z, 32) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", Z[3], Z[2], Z[1], Z[0]);
            fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", z[3], z[2], z[1], z[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_subp(void) {
    FILE* file = fopen("tests/subp.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], z[4], Z[4];

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", z+3, z+2, z+1, z+0); fgetc(file);

        prime_sub(Z, x, y);

        if (memcmp(Z, z, 32) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", Z[3], Z[2], Z[1], Z[0]);
            fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", z[3], z[2], z[1], z[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_mulp(void) {
    FILE* file = fopen("tests/mulp.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], z[4], Z[4];

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", z+3, z+2, z+1, z+0); fgetc(file);

        prime_mul(Z, x, y);
        prime_mod(Z); /* Add final mod operation, because i haven't results for reduced mul */

        if (memcmp(Z, z, 32) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", Z[3], Z[2], Z[1], Z[0]);
            fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", z[3], z[2], z[1], z[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_addo(void) {
    FILE* file = fopen("tests/addo.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], z[4], Z[4];

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", z+3, z+2, z+1, z+0); fgetc(file);

        order_add(Z, x, y);

        if (memcmp(Z, z, 32) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", Z[3], Z[2], Z[1], Z[0]);
            fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", z[3], z[2], z[1], z[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_mulo(void) {
    FILE* file = fopen("tests/mulo.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], z[4], Z[4];

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", z+3, z+2, z+1, z+0); fgetc(file);

        order_mul(Z, x, y);
        order_mod(Z); /* Add final mod operation, because i haven't results for reduced mul */

        if (memcmp(Z, z, 32) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", Z[3], Z[2], Z[1], Z[0]);
            fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", z[3], z[2], z[1], z[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_exp(void) {
    FILE* file = fopen("tests/exp.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], Y[4];

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);

        prime_exp_by_2_252_3(Y, x);
        prime_mod(Y); /* Add final mod operation, because i haven't results for reduced mul */

        if (memcmp(Y, y, 32) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", Y[3], Y[2], Y[1], Y[0]);
            fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", y[3], y[2], y[1], y[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_inv(void) {
    FILE* file = fopen("tests/inv.txt", "r"); size_t i, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint64_t x[4], y[4], Y[4];

        fscanf(file, "%16lx%16lx%16lx%16lx", x+3, x+2, x+1, x+0); fgetc(file);
        fscanf(file, "%16lx%16lx%16lx%16lx", y+3, y+2, y+1, y+0); fgetc(file);

        prime_inv(Y, x);
        prime_mod(Y); /* Add final mod operation, because i haven't results for reduced mul */

        if (memcmp(Y, y, 32) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", Y[3], Y[2], Y[1], Y[0]);
            fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", y[3], y[2], y[1], y[0]);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_pack(void) {
    static const uint64_t gp[4] = { 0x6666666666666658, 0x6666666666666666, 0x6666666666666666, 0x6666666666666666 };
    static const uint64_t g[2][4] = {
        { 0xc9562d608f25d51a, 0x692cc7609525a7b2, 0xc0a4e231fdd6dc5c, 0x216936d3cd6e53fe },
        { 0x6666666666666658, 0x6666666666666666, 0x6666666666666666, 0x6666666666666666 }
    };
    uint64_t Gp[4], G[2][4]; int e = 0;

    point_pack(Gp, g);
    point_unpack(G, gp);

    if (memcmp(Gp, gp, 32) != 0) {
        fprintf(stderr, "Packing\n");
        fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", Gp[3], Gp[2], Gp[1], Gp[0]);
        fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", gp[3], gp[2], gp[1], gp[0]);
        ++e;
    }
    if (memcmp(G[0], g[0], 32) != 0) {
        fprintf(stderr, "Unpacking x\n");
        fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", G[0][3], G[0][2], G[0][1], G[0][0]);
        fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", g[0][3], g[0][2], g[0][1], g[0][0]);
        ++e;
    }
    if (memcmp(G[1], g[1], 32) != 0) {
        fprintf(stderr, "Unpacking y\n");
        fprintf(stderr, "  wrong  : %016lx%016lx%016lx%016lx\n", G[1][3], G[1][2], G[1][1], G[1][0]);
        fprintf(stderr, "  correct: %016lx%016lx%016lx%016lx\n", g[1][3], g[1][2], g[1][1], g[1][0]);
        ++e;
    }

    return e > 0;
}

int test_sha(void) {
    FILE* file = fopen("tests/sha.txt", "r"); size_t i, j, e;
    if (!file) return 1;

    for (i = e = 0; i < 100; i++) {
        uint8_t c[209], m[99], h[64], H[64];

        for (j = 0; j <  i; j++) { m[j] = fgetc(file); } fgetc(file);
        for (j = 0; j < 64; j++) { unsigned b; fscanf(file, "%2x", &b); h[j] = b; } fgetc(file);

        shahf_begin(c);
        shahf_update(c, m, i);
        shahf_end(c, H);

        if (memcmp(H, h, 64) != 0) {
            fprintf(stderr, "Case #%lu\n", i);
            fprintf(stderr, "  wrong  : "); for (j = 0; j < 64; j++) { fprintf(stderr, "%02x", H[j]); } fputc('\n', stderr);
            fprintf(stderr, "  correct: "); for (j = 0; j < 64; j++) { fprintf(stderr, "%02x", h[j]); } fputc('\n', stderr);
            if (++e >= 10) break;
        }
    }

    fclose(file);
    return e > 0;
}

int test_keygen(void) {
    const uint8_t pvtk[32] = {
        "\x9d\x61\xb1\x9d\xef\xfd\x5a\x60\xba\x84\x4a\xf4\x92\xec\x2c\xc4"
        "\x44\x49\xc5\x69\x7b\x32\x69\x19\x70\x3b\xac\x03\x1c\xae\x7f\x60" };
    const uint8_t pubk[32] = {
        "\xd7\x5a\x98\x01\x82\xb1\x0a\xb7\xd5\x4b\xfe\xd3\xc9\x64\x07\x3a"
        "\x0e\xe1\x72\xf3\xda\xa6\x23\x25\xaf\x02\x1a\x68\xf7\x07\x51\x1a" };
    uint8_t PvtK[32], PubK[32]; size_t j, e = 0;

    eddsa_ed25519_key_generation(PvtK, PubK, pvtk);

    if (memcmp(PvtK, pvtk, 32) != 0) {
        fprintf(stderr, "Private key\n");
        fprintf(stderr, "  wrong  : "); for (j = 0; j < 32; j++) { fprintf(stderr, "%02x", PvtK[j]); } fputc('\n', stderr);
        fprintf(stderr, "  correct: "); for (j = 0; j < 32; j++) { fprintf(stderr, "%02x", pvtk[j]); } fputc('\n', stderr);
        ++e;
    }
    if (memcmp(PubK, pubk, 32) != 0) {
        fprintf(stderr, "Public key\n");
        fprintf(stderr, "  wrong  : "); for (j = 0; j < 32; j++) { fprintf(stderr, "%02x", PubK[j]); } fputc('\n', stderr);
        fprintf(stderr, "  correct: "); for (j = 0; j < 32; j++) { fprintf(stderr, "%02x", pubk[j]); } fputc('\n', stderr);
        ++e;
    }

    return e > 0;
}

#define test_case(name) do { \
    int rc = test_ ## name(); \
    printf("Run test " #name " ... " \
        "\x1b[3%s\x1b[m\n", rc ? "1mFAIL" : "2mOK"); \
    if (rc) return 1; \
} while (0)

int main(void) {
    test_case(add);
    test_case(sub);
    test_case(mul);

    test_case(addp);
    test_case(subp);
    test_case(mulp);

    test_case(addo);
    test_case(mulo);

    test_case(exp);
    test_case(inv);

    test_case(pack);

    test_case(sha);

    test_case(keygen);
    return 0;
}
