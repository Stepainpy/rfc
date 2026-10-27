#include "eddsa_ed25519.h"

/* Adaptation of code from GitHub repo msotoodeh/curve25519 */

#include <string.h>

/* * * * * * * * * * * * * * * * * *
 * Basic multiprecision arithmetic *
 * * * * * * * * * * * * * * * * * */

#  if defined(__GNUC__)

#include <stdint.h>

typedef uint8_t u8_t;
typedef uint64_t u64_t;

#define instr_swp __builtin_bswap64

u64_t instr_ror(u64_t n, int s) {
    s &= 63; return n >> s | n << (-s & 63);
}

#if ~0ul == UINT64_MAX
#  define instr_adc __builtin_addcl
#  define instr_sbb __builtin_subcl
#else
#  define instr_adc __builtin_addcll
#  define instr_sbb __builtin_subcll
#endif

u64_t instr_mul(u64_t x, u64_t y, u64_t* h) {
    __uint128_t z = (__uint128_t)x * (__uint128_t)y;
    *h = z >> 64; return z;
}

#elif defined(_MSC_VER)

#include <stdlib.h>
#include <intrin.h>

typedef unsigned char u8_t;
typedef unsigned __int64 u64_t;

#define instr_swp _byteswap_uint64
#define instr_ror _rotr64

u64_t instr_adc(u64_t x, u64_t y, u64_t ci, u64_t* co) {
    u64_t z; *co = _addcarry_u64(ci, x, y, &z);
    return z;
}

u64_t instr_sbb(u64_t x, u64_t y, u64_t ci, u64_t* co) {
    u64_t z; *co = _subborrow_u64(ci, x, y, &z);
    return z;
}

#define instr_mul _umul128

#else /* neither GNU nor Clang nor MSVC */
#  error "Unsupport compiler"
#endif

u64_t arith_add(u64_t z[4], const u64_t x[4], const u64_t y[4]) {
    u64_t c = 0;
    z[0] = instr_adc(x[0], y[0], c, &c);
    z[1] = instr_adc(x[1], y[1], c, &c);
    z[2] = instr_adc(x[2], y[2], c, &c);
    z[3] = instr_adc(x[3], y[3], c, &c);
    return c; /* 0 or 1 */
}

u64_t arith_sub(u64_t z[4], const u64_t x[4], const u64_t y[4]) {
    u64_t b = 0;
    z[0] = instr_sbb(x[0], y[0], b, &b);
    z[1] = instr_sbb(x[1], y[1], b, &b);
    z[2] = instr_sbb(x[2], y[2], b, &b);
    z[3] = instr_sbb(x[3], y[3], b, &b);
    return -b; /* 0 or ~0 */
}

void prime_add(u64_t z[4], const u64_t x[4], const u64_t y[4]) {
    u64_t t[4] = {0};

    t[0] = 38 * arith_add(z, x, y);
    t[0] = 38 * arith_add(z, z, t);
                arith_add(z, z, t);
}

void prime_sub(u64_t z[4], const u64_t x[4], const u64_t y[4]) {
    u64_t t[4] = {0};

    t[0] = 38 & arith_sub(z, x, y);
    t[0] = 38 & arith_sub(z, z, t);
                arith_sub(z, z, t);
}

void arith_mul_set(u64_t y[5], u64_t n, const u64_t x[4]) {
    u64_t l, h, c, hc;

    l = instr_mul(n, x[0], &h);
    y[0] = l; hc = h;

    l = instr_mul(n, x[1], &h);
    l = instr_adc(l, hc, 0, &c); h += c;
    y[1] = l; hc = h;

    l = instr_mul(n, x[2], &h);
    l = instr_adc(l, hc, 0, &c); h += c;
    y[2] = l; hc = h;

    l = instr_mul(n, x[3], &h);
    l = instr_adc(l, hc, 0, &c); h += c;
    y[3] = l; hc = h;

    y[4] = hc;
}

void arith_mul_add(u64_t y[5], u64_t n, const u64_t x[4]) {
    u64_t l, h, c, hc;

    l = instr_mul(n, x[0], &h);
    l = instr_adc(l, y[0], 0, &c); h += c;
    y[0] = l; hc = h;

    l = instr_mul(n, x[1], &h);
    l = instr_adc(l, y[1], 0, &c); h += c;
    l = instr_adc(l,  hc , 0, &c); h += c;
    y[1] = l; hc = h;

    l = instr_mul(n, x[2], &h);
    l = instr_adc(l, y[2], 0, &c); h += c;
    l = instr_adc(l,  hc , 0, &c); h += c;
    y[2] = l; hc = h;

    l = instr_mul(n, x[3], &h);
    l = instr_adc(l, y[3], 0, &c); h += c;
    l = instr_adc(l,  hc , 0, &c); h += c;
    y[3] = l; hc = h;

    y[4] = hc;
}

void arith_mul(u64_t z[8], const u64_t x[4], const u64_t y[4]) {
    arith_mul_set(z + 0, x[0], y);
    arith_mul_add(z + 1, x[1], y);
    arith_mul_add(z + 2, x[2], y);
    arith_mul_add(z + 3, x[3], y);
}

void arith_wma_rdc(u64_t z[4], u64_t n, const u64_t x[4], const u64_t y[4]) {
    u64_t l, h, c, hc;

    l = instr_mul(n, x[0], &h);
    l = instr_adc(l, y[0], 0, &c); h += c;
    z[0] = l; hc = h;

    l = instr_mul(n, x[1], &h);
    l = instr_adc(l, y[1], 0, &c); h += c;
    l = instr_adc(l,  hc , 0, &c); h += c;
    z[1] = l; hc = h;

    l = instr_mul(n, x[2], &h);
    l = instr_adc(l, y[2], 0, &c); h += c;
    l = instr_adc(l,  hc , 0, &c); h += c;
    z[2] = l; hc = h;

    l = instr_mul(n, x[3], &h);
    l = instr_adc(l, y[3], 0, &c); h += c;
    l = instr_adc(l,  hc , 0, &c); h += c;
    z[3] = l; hc = h;

    l = instr_mul(38, hc, &h);
    z[0] = instr_adc(z[0], l, 0, &c);
    z[1] = instr_adc(z[1], h, c, &c);
    z[2] = instr_adc(z[2], 0, c, &c);
    z[3] = instr_adc(z[3], 0, c, &c);

    l = 38 * c;
    z[0] = instr_adc(z[0], l, c, &c);
    z[1] = instr_adc(z[1], 0, c, &c);
    z[2] = instr_adc(z[2], 0, c, &c);
    z[3] = instr_adc(z[3], 0, c, &c);
}

void prime_mul(u64_t z[4], const u64_t x[4], const u64_t y[4]) {
    u64_t t[8];
    arith_mul(t, x, y);
    arith_wma_rdc(z, 38, t + 4, t);
}

void prime_sqr(u64_t y[4], const u64_t x[4]) {
    u64_t t[8];
    arith_mul(t, x, x);
    arith_wma_rdc(y, 38, t + 4, t);
}

static const u64_t const_prime[4] = {
    0xffffffffffffffed, 0xffffffffffffffff,
    0xffffffffffffffff, 0x7fffffffffffffff
};

void prime_mod(u64_t x[4]) {
    u64_t t[4], m;

    m = arith_sub(x, x, const_prime);
    t[0] = m & ~(u64_t)18; t[1] = t[2] = m; t[3] = m >> 1;
    arith_add(x, x, t);

    m = arith_sub(x, x, const_prime);
    t[0] = m & ~(u64_t)18; t[1] = t[2] = m; t[3] = m >> 1;
    arith_add(x, x, t);
}

static const u64_t const_minus_r[4] = {
    0x812631a5cf5d3ed0, 0x4def9dea2f79cd65, 1, 0
};

static const u64_t const_nxo[16][4] = {
    { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000 },
    { 0x5812631a5cf5d3ed, 0x14def9dea2f79cd6, 0x0000000000000000, 0x1000000000000000 },
    { 0xb024c634b9eba7da, 0x29bdf3bd45ef39ac, 0x0000000000000000, 0x2000000000000000 },
    { 0x0837294f16e17bc7, 0x3e9ced9be8e6d683, 0x0000000000000000, 0x3000000000000000 },
    { 0x60498c6973d74fb4, 0x537be77a8bde7359, 0x0000000000000000, 0x4000000000000000 },
    { 0xb85bef83d0cd23a1, 0x685ae1592ed6102f, 0x0000000000000000, 0x5000000000000000 },
    { 0x106e529e2dc2f78e, 0x7d39db37d1cdad06, 0x0000000000000000, 0x6000000000000000 },
    { 0x6880b5b88ab8cb7b, 0x9218d51674c549dc, 0x0000000000000000, 0x7000000000000000 },
    { 0xc09318d2e7ae9f68, 0xa6f7cef517bce6b2, 0x0000000000000000, 0x8000000000000000 },
    { 0x18a57bed44a47355, 0xbbd6c8d3bab48389, 0x0000000000000000, 0x9000000000000000 },
    { 0x70b7df07a19a4742, 0xd0b5c2b25dac205f, 0x0000000000000000, 0xa000000000000000 },
    { 0xc8ca4221fe901b2f, 0xe594bc9100a3bd35, 0x0000000000000000, 0xb000000000000000 },
    { 0x20dca53c5b85ef1c, 0xfa73b66fa39b5a0c, 0x0000000000000000, 0xc000000000000000 },
    { 0x78ef0856b87bc309, 0x0f52b04e4692f6e2, 0x0000000000000001, 0xd000000000000000 },
    { 0xd1016b71157196f6, 0x2431aa2ce98a93b8, 0x0000000000000001, 0xe000000000000000 },
    { 0x2913ce8b72676ae3, 0x3910a40b8c82308f, 0x0000000000000001, 0xf000000000000000 },
};

void arith_hw_rdc(u64_t y[4], u64_t n, const u64_t x[4]) {
    u64_t t[4], h, l, hc, c;

    l = instr_mul(n, const_minus_r[0], &h);
    t[0] = l; hc = h;

    l = instr_mul(n, const_minus_r[1], &h);
    l = instr_adc(l, hc, 0, &c); h += c;
    t[1] = l; hc = h;

    l = instr_mul(n, const_minus_r[2], &h);
    l = instr_adc(l, hc, 0, &c); h += c;
    t[2] = l;
    t[3] = h;

    c = arith_sub(y, x, t);
    arith_add(y, y, const_nxo[c & 1]);
}

void order_from_digest(u64_t x[4], const void* digest) {
    u64_t t[8]; memcpy(t, digest, 64);
    arith_hw_rdc(t + 3, t[7], t + 3);
    arith_hw_rdc(t + 2, t[6], t + 2);
    arith_hw_rdc(t + 1, t[5], t + 1);
    arith_hw_rdc(  x  , t[4], t + 0);
}

void order_add(u64_t z[4], const u64_t x[4], const u64_t y[4]) {
    u64_t c = arith_add(z, x, y);
    arith_hw_rdc(z, c, z);
}

void order_mul(u64_t z[4], const u64_t x[4], const u64_t y[4]) {
    u64_t t[8]; arith_mul(t, x, y);
    order_from_digest(z, t);
}

void order_mod(u64_t x[4]) {
    u64_t c = arith_sub(x, x, const_nxo[x[3] >> 60]);
    arith_add(x, x, const_nxo[c & 1]);
}

/* * * * * * * * * * * * * * * * * * *
 * Complex multiprecision arithmetic *
 * * * * * * * * * * * * * * * * * * */

void prime_inv(u64_t y[4], const u64_t x[4]) {
    u64_t x2_5_0[4], x2_10_0[4], x2_20_0[4], x2_50_0[4], x2_100_0[4];
    u64_t t0[4], t1[4], x2[4], x9[4], x11[4], i;

    prime_sqr(x2, x);
    prime_sqr(t1, x2);
    prime_sqr(t0, t1);
    prime_mul(x9, t0, x);
    prime_mul(x11, x9, x2);
    prime_sqr(t0, x11);
    prime_mul(x2_5_0, t0, x9);

    prime_sqr(t0, x2_5_0);
    prime_sqr(t1, t0);
    prime_sqr(t0, t1);
    prime_sqr(t1, t0);
    prime_sqr(t0, t1);
    prime_mul(x2_10_0, t0, x2_5_0);

    prime_sqr(t0, x2_10_0);
    prime_sqr(t1, t0);
    for (i = 2; i < 10; i += 2)
        prime_sqr(t0, t1),
        prime_sqr(t1, t0);
    prime_mul(x2_20_0, t1, x2_10_0);

    prime_sqr(t0, x2_20_0);
    prime_sqr(t1, t0);
    for (i = 2; i < 20; i += 2)
        prime_sqr(t0, t1),
        prime_sqr(t1, t0);
    prime_mul(t0, t1, x2_20_0);

    prime_sqr(t1, t0);
    prime_sqr(t0, t1);
    for (i = 2; i < 10; i += 2)
        prime_sqr(t1, t0),
        prime_sqr(t0, t1);
    prime_mul(x2_50_0, t0, x2_10_0);

    prime_sqr(t0, x2_50_0);
    prime_sqr(t1, t0);
    for (i = 2; i < 50; i += 2)
        prime_sqr(t0, t1),
        prime_sqr(t1, t0);
    prime_mul(x2_100_0, t1, x2_50_0);

    prime_sqr(t1, x2_100_0);
    prime_sqr(t0, t1);
    for (i = 2; i < 100; i += 2)
        prime_sqr(t1, t0),
        prime_sqr(t0, t1);
    prime_mul(t1, t0, x2_100_0);

    prime_sqr(t0, t1);
    prime_sqr(t1, t0);
    for (i = 2; i < 50; i += 2)
        prime_sqr(t0, t1),
        prime_sqr(t1, t0);
    prime_mul(t0, t1, x2_50_0);

    prime_sqr(t1, t0);
    prime_sqr(t0, t1);
    prime_sqr(t1, t0);
    prime_sqr(t0, t1);
    prime_sqr(t1, t0);

    prime_mul(y, t1, x11);
}

void prime_sqr_by_n(u64_t z[4], const u64_t x[4], u64_t n, const u64_t y[4]) {
    u64_t t[4];
    prime_sqr(t, x);
    while (n --> 1) prime_sqr(t, t);
    prime_mul(z, t, y);
}

void prime_exp_by_2_252_3(u64_t y[4], const u64_t x[4]) {
    u64_t x2[4], x5[4], x9[4], x10[4], x11[4];
    u64_t x20[4], x50[4], x100[4], t[4];

    prime_sqr(x2, x);
    prime_sqr_by_n(x9, x2, 2, x);
    prime_mul(x11, x9, x2);
    prime_sqr(t, x11);
    prime_mul(x5, t, x9);
    prime_sqr_by_n(x10, x5, 5, x5);
    prime_sqr_by_n(x20, x10, 10, x10);
    prime_sqr_by_n(t, x20, 20, x20);
    prime_sqr_by_n(x50, t, 10, x10);
    prime_sqr_by_n(x100, x50, 50, x50);
    prime_sqr_by_n(t, x100, 100, x100);
    prime_sqr_by_n(t, t, 50, x50);
    prime_sqr(t, t);
    prime_sqr(t, t);
    prime_mul(y, t, x);
}

/* * * * * * * * * * * *
 * Operation on points *
 * * * * * * * * * * * */

typedef struct { u64_t x[4], y[4]; } point_affine_t;
typedef struct { u64_t x[4], y[4], z[4], t[4]; } point_proj_t;

typedef struct {
    u64_t ypx[4]; /* y + x */
    u64_t ymx[4]; /* y - x */
    u64_t t2d[4]; /* 2d * t */
} point_pc_affine_t;

typedef struct {
    u64_t ypx[4]; /* y + x */
    u64_t ymx[4]; /* y - x */
    u64_t t2d[4]; /* 2d * t */
    u64_t z2 [4]; /* 2 * z */
} point_pc_proj_t;

void point_pack(void* out, const point_affine_t* p) {
    u64_t* o = out;

    o[0] = p->y[0]; o[1] = p->y[1]; o[2] = p->y[2];
    o[3] = (p->y[3] & (~(u64_t)0 >> 1)) | p->x[0] << 63;
}

static const u64_t const_d[4] = {
    0x75eb4dca135978a3, 0x00700a4d4141d8ab,
    0x8cc740797779e898, 0x52036cee2b6ffe73
};

static const u64_t const_i[4] = {
    0xc4ee1b274a0ea0b0, 0x2f431806ad2fe478,
    0x2b4d00993dfbd7a7, 0x2b8324804fc1df0b
};

void point_unpack(point_affine_t* p, const void* in) {
    static const u64_t zero[4] = { 0, 0, 0, 0 };
    static const u64_t one [4] = { 1, 0, 0, 0 };

    u64_t u[4], v[4], a[4], t[4], parity;

    memcpy(p->y, in, 32);
    parity = p->y[3] >> 63;
    p->y[3] &= ~(u64_t)0 >> 1;

    prime_sqr(u, p->y);
    prime_mul(v, u, const_d);
    prime_sub(u, u, one);
    prime_add(v, v, one);

    /* use more short equation from errata-5758 */
    prime_mul(a, u, v);
    prime_exp_by_2_252_3(a, a);
    prime_mul(p->x, u, a);

    prime_sqr(a, p->x);
    prime_mul(a, a, v);
    prime_sub(a, a, u);
    prime_mod(a);
    if (memcmp(a, zero, 32) != 0)
        prime_mul(p->x, p->x, const_i);

    while (arith_sub(t, p->x, const_prime) == 0)
        memcpy(p->x, t, 32);

    if ((p->x[0] & 1) != parity)
        arith_sub(p->x, const_prime, p->x);
}

void point_add(point_proj_t* p, const point_pc_affine_t* q) {
    u64_t a[4], b[4], c[4], d[4], e[4];

    prime_sub(a, p->y, p->x);
    prime_mul(a, a, q->ymx);
    prime_add(b, p->y, p->x);
    prime_mul(b, b, q->ypx);
    prime_mul(c, p->t, q->t2d);
    prime_add(d, p->z, p->z);

    prime_sub(e, b, a);
    prime_add(b, b, a);
    prime_sub(a, d, c);
    prime_add(d, d, c);

    prime_mul(p->x, e, a);
    prime_mul(p->y, b, d);
    prime_mul(p->z, d, a);
    prime_mul(p->t, e, b);
}

static const u64_t const_double_prime[4] = {
    0xffffffffffffffda, 0xffffffffffffffff,
    0xffffffffffffffff, 0xffffffffffffffff
};

void point_dbl(point_proj_t* p) {
    u64_t a[4], b[4], c[4], d[4], e[4];

    prime_sqr(a, p->x);
    prime_sqr(b, p->y);
    prime_sqr(c, p->z);
    prime_add(c, c, c);
    prime_sub(d, const_double_prime, a);

    prime_sub(a, d, b);
    prime_add(d, d, b);
    prime_sub(b, d, c);
    prime_add(e, p->x, p->y);
    prime_sqr(e, e);
    prime_add(e, e, a);

    prime_mul(p->x, e, b);
    prime_mul(p->y, a, d);
    prime_mul(p->z, d, b);
    prime_mul(p->t, e, a);
}

/* * * * * * * * * * * * * * * *
 * SHA-2 512-bit hash function *
 * * * * * * * * * * * * * * * */

typedef struct {
    u64_t h[8], len[2];
    u8_t buf[128], cnt;
} shahf_ctx_t;

static const u64_t const_k[80] = {
    0x428a2f98d728ae22, 0x7137449123ef65cd, 0xb5c0fbcfec4d3b2f, 0xe9b5dba58189dbbc, 0x3956c25bf348b538,
    0x59f111f1b605d019, 0x923f82a4af194f9b, 0xab1c5ed5da6d8118, 0xd807aa98a3030242, 0x12835b0145706fbe,
    0x243185be4ee4b28c, 0x550c7dc3d5ffb4e2, 0x72be5d74f27b896f, 0x80deb1fe3b1696b1, 0x9bdc06a725c71235,
    0xc19bf174cf692694, 0xe49b69c19ef14ad2, 0xefbe4786384f25e3, 0x0fc19dc68b8cd5b5, 0x240ca1cc77ac9c65,
    0x2de92c6f592b0275, 0x4a7484aa6ea6e483, 0x5cb0a9dcbd41fbd4, 0x76f988da831153b5, 0x983e5152ee66dfab,
    0xa831c66d2db43210, 0xb00327c898fb213f, 0xbf597fc7beef0ee4, 0xc6e00bf33da88fc2, 0xd5a79147930aa725,
    0x06ca6351e003826f, 0x142929670a0e6e70, 0x27b70a8546d22ffc, 0x2e1b21385c26c926, 0x4d2c6dfc5ac42aed,
    0x53380d139d95b3df, 0x650a73548baf63de, 0x766a0abb3c77b2a8, 0x81c2c92e47edaee6, 0x92722c851482353b,
    0xa2bfe8a14cf10364, 0xa81a664bbc423001, 0xc24b8b70d0f89791, 0xc76c51a30654be30, 0xd192e819d6ef5218,
    0xd69906245565a910, 0xf40e35855771202a, 0x106aa07032bbd1b8, 0x19a4c116b8d2d0c8, 0x1e376c085141ab53,
    0x2748774cdf8eeb99, 0x34b0bcb5e19b48a8, 0x391c0cb3c5c95a63, 0x4ed8aa4ae3418acb, 0x5b9cca4f7763e373,
    0x682e6ff3d6b2b8a3, 0x748f82ee5defb2fc, 0x78a5636f43172f60, 0x84c87814a1f0ab72, 0x8cc702081a6439ec,
    0x90befffa23631e28, 0xa4506cebde82bde9, 0xbef9a3f7b2c67915, 0xc67178f2e372532b, 0xca273eceea26619c,
    0xd186b8c721c0c207, 0xeada7dd6cde0eb1e, 0xf57d4f7fee6ed178, 0x06f067aa72176fba, 0x0a637dc5a2c898a6,
    0x113f9804bef90dae, 0x1b710b35131c471b, 0x28db77f523047d84, 0x32caab7b40c72493, 0x3c9ebe0a15c9bebc,
    0x431d67c49c100d4c, 0x4cc5d4becb3e42b6, 0x597f299cfc657e2a, 0x5fcb6fab3ad6faec, 0x6c44198c4a475817,
};

void shahf_round(shahf_ctx_t* ctx) {
    u64_t a, b, c, d, e, f, g, h;
    u64_t s0, s1, ch, ma, t1, t2;
    u64_t w[80]; size_t i;

    memcpy(w, ctx->buf, 128);
    memset(ctx->buf, 0, 128);
    ctx->cnt = 0;

    for (i = 0; i < 16; i++)
        w[i] = instr_swp(w[i]);
    for (i = 16; i < 80; i++) {
        s0 = instr_ror(w[i - 15],  1) ^ instr_ror(w[i - 15],  8) ^ (w[i - 15] >> 7);
        s1 = instr_ror(w[i -  2], 19) ^ instr_ror(w[i -  2], 61) ^ (w[i -  2] >> 6);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }

    a = ctx->h[0]; b = ctx->h[1]; c = ctx->h[2]; d = ctx->h[3];
    e = ctx->h[4]; f = ctx->h[5]; g = ctx->h[6]; h = ctx->h[7];

    for (i = 0; i < 80; i++) {
        s0 = instr_ror(a, 28) ^ instr_ror(a, 34) ^ instr_ror(a, 39);
        s1 = instr_ror(e, 14) ^ instr_ror(e, 18) ^ instr_ror(e, 41);
        ch = (e & f) ^ (~e & g);
        ma = (a & b) ^ (a & c) ^ (b & c);
        t1 = h + s1 + ch + const_k[i] + w[i];
        t2 = s0 + ma;

        h = g; g = f; f = e; e = t1 + d;
        d = c; c = b; b = a; a = t1 + t2;
    }

    ctx->h[0] += a; ctx->h[1] += b; ctx->h[2] += c; ctx->h[3] += d;
    ctx->h[4] += e; ctx->h[5] += f; ctx->h[6] += g; ctx->h[7] += h;
}

void shahf_begin(shahf_ctx_t* ctx) {
    memset(ctx, 0, sizeof *ctx);
    ctx->h[0] = 0x6a09e667f3bcc908;
    ctx->h[1] = 0xbb67ae8584caa73b;
    ctx->h[2] = 0x3c6ef372fe94f82b;
    ctx->h[3] = 0xa54ff53a5f1d36f1;
    ctx->h[4] = 0x510e527fade682d1;
    ctx->h[5] = 0x9b05688c2b3e6c1f;
    ctx->h[6] = 0x1f83d9abfb41bd6b;
    ctx->h[7] = 0x5be0cd19137e2179;
}

void shahf_update(shahf_ctx_t* ctx, const void* data, size_t size) {
    while (size) {
        size_t rem = 128 - ctx->cnt; u64_t c;
        size_t min = size < rem ? size : rem;

        memcpy(ctx->buf + ctx->cnt, data, min);
        data = (const u8_t*)data + min;
        ctx->cnt += min; size -= min;

        ctx->len[0] = instr_adc(ctx->len[0], 8 * min, 0, &c);
        ctx->len[1] = instr_adc(ctx->len[1],    0   , c, &c);

        if (ctx->cnt == 128)
            shahf_round(ctx);
    }
}

void shahf_end(shahf_ctx_t* ctx, void* hash) {
    size_t i; ctx->buf[ctx->cnt++] = 0x80;
    if (ctx->cnt > 112) shahf_round(ctx);

    ((u64_t*)ctx->buf)[14] = instr_swp(ctx->len[1]);
    ((u64_t*)ctx->buf)[15] = instr_swp(ctx->len[0]);
    shahf_round(ctx);

    for (i = 0; i < 8; i++)
        ctx->h[i] = instr_swp(ctx->h[i]);
    memcpy(hash, ctx->h, 64);
}

/* * * * * * * * * * *
 * Ed25519 functions *
 * * * * * * * * * * */

void curve_fold8(u8_t f[32], const u64_t x[4]) {
    u64_t t[4]; size_t i, j;
    memcpy(t, x, 32);
    for (i = 0; i < 32; i++) {
        f[i] = 0;
        for (j = 4; j --> 0;) {
            f[i] = f[i] << 1 | (t[j] >> 63 & 1);
            f[i] = f[i] << 1 | (t[j] >> 31 & 1);
            t[j] <<= 1;
        }
    }
}

static const point_pc_affine_t const_base_fold[256] = {
       {
        { 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, },
        { 0x0000000000000001, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, },
        { 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, 0x0000000000000000, },
    }, {
        { 0x2fbc93c6f58c3b85, 0xcf932dc6fb8c0e19, 0x270b4898643d42c2, 0x07cf9d3a33d4ba65, },
        { 0x9d103905d740913e, 0xfd399f05d140beb3, 0xa5c18434688f8a09, 0x44fd2f9298f81267, },
        { 0xabc91205877aaa68, 0x26d9e823ccaac49e, 0x5a1b7dcbdd43598c, 0x6f117b689f0c65a8, },
    }, {
        { 0x8765b69f7b85c5e8, 0x6ff0678bd168bab2, 0x3a70e77c1d330f9b, 0x3a5f6d51b0af8e7c, },
        { 0x61368756a60dac5f, 0x17e02f6aebabdc57, 0x7f193f2d4cce0f7d, 0x20234a7789ecdcf0, },
        { 0x76d20db67178b252, 0x071c34f9d51ed160, 0xf62a4a20b3e41170, 0x7cd682353cffe366, },
    }, {
        { 0x812d6a74a76c971f, 0x1ceca25f0739f84d, 0x27d06c3c4102b810, 0x7edc1cc6031906f7, },
        { 0x500519ad43ba112a, 0xebda0d43aebd3ad9, 0xf3f6b90c3e2cf3f0, 0x0de1829ba60cd4f5, },
        { 0xf5472c0e2ab6034b, 0xbd4dbf55e480f342, 0xd443d2710cced266, 0x4fe5ca429c2db4f7, },
    }, {
        { 0xcd2a65e777d1f515, 0x548991878faa60f1, 0xb1b73bbcdabc06e5, 0x654878cba97cc9fb, },
        { 0x51138ec78df6b0fe, 0x5397da89e575f51b, 0x09207a1d717af1b9, 0x2102fdba2b20d650, },
        { 0x969ee405055ce6a1, 0x36bca7681251ad29, 0x3a1af517aa7da415, 0x0ad725db29ecb2ba, },
    }, {
        { 0x0055c585601e59e8, 0x8793342b66480e60, 0x3e14aad0fe45e44c, 0x26ead8e64813cf2b, },
        { 0xcb75b8b69c8462a4, 0x2dd86fc567d31cd7, 0xcd1972ec881342f6, 0x0975b5970fc12f2f, },
        { 0x63cf2303da5ba743, 0x04bf9d8152f1ba6e, 0x333790d0aa7367da, 0x534670479df6c5ea, },
    }, {
        { 0x0af0993c7738657c, 0xa9fb4a4747181bd2, 0x55dde7ab3e9ec14a, 0x0b9ec7f637c8321c, },
        { 0x0ae97f7d88035653, 0x17fb7b108d25068e, 0x6328fc0e3c1e8c9a, 0x5441585e50fe219b, },
        { 0xcb52864c640f0146, 0x94562266317ef0d3, 0xb876f4c41595bb5d, 0x4ed2f874af7f7713, },
    }, {
        { 0xfeb687d874d2cda7, 0xb8042dc1ae27fb22, 0x1564556d2c93adcc, 0x11442d28e3a17816, },
        { 0xcec3358507f0c5d4, 0xfb9e7ef995b01c20, 0x0a869383dea13143, 0x16c63a61f5f8c59b, },
        { 0xded6638cd6fed546, 0x1abdb60b3803190e, 0xe0c56f156a614da3, 0x04c548569edc6e88, },
    }, {
        { 0xd598639c12ddb0a4, 0xa5d19f30c024866b, 0xd17c2f0358fce460, 0x07a195152e095e8a, },
        { 0x296fa9c59c2ec4de, 0xbc8b61bf4f84f3cb, 0x1c7706d917a8f908, 0x63b795fc7ad3255d, },
        { 0xa8368f02389e5fc8, 0x90433b02cf8de43b, 0xafa1fd5dc5412643, 0x3e8fe83d032f0137, },
    }, {
        { 0xb5f7f69a17d3a339, 0x805e614536a01f1e, 0x08e62b6d9b01a221, 0x0af8325038f1898f, },
        { 0xbcc2be89bdf01e71, 0x6758f9e1c3280d0f, 0xf2a56cdeccf36c58, 0x1cdfbf7e80faacc0, },
        { 0x730688fb39bfb3b2, 0xf493c37607c06f81, 0xbda608c28db1ba83, 0x2990e0a263ecca60, },
    }, {
        { 0xfca8ea4159dc4791, 0x0fae3dab8b3aa058, 0xbe13396f4ee996eb, 0x379d09bb51936c6f, },
        { 0xb6601a1bf614affb, 0x8360c886210392ea, 0x4867333c56349198, 0x03224a6ff049c42c, },
        { 0x5e267a306fb88974, 0xbeb84f824f8fb990, 0x6029b7b918c57b0d, 0x60670bbea2357df1, },
    }, {
        { 0xc8e97558ec02dcd8, 0x70e6d52d734215b0, 0x5ad04683c97f37cc, 0x2c610ab3c1d685ab, },
        { 0xabc333f4e10c024f, 0x4f8149e4672386cc, 0x3bad171c4e7bf204, 0x4cea190be78f5067, },
        { 0x8d8d3a7902f3325f, 0x549cbd48bf8abde9, 0x9da1a0f611d19a0c, 0x718bf518372c817c, },
    }, {
        { 0x21c5e9f32381ad37, 0xc340283eb96cd050, 0xfbdee4713150abca, 0x31b24deec8fc9db7, },
        { 0x2c160eed2b76ccbb, 0x6dd5c36b54d82775, 0x7f278210d5e0a5e9, 0x076fc47cc59cf540, },
        { 0xbdf661d85f61df81, 0xfbb57fe54eaec5bd, 0x5aa7105ee97cf86e, 0x00a448b6d8d3cd95, },
    }, {
        { 0x584444a6972952b6, 0x1a2e4af7f3ba527c, 0x5ad4fff4c1e4ed75, 0x4c05a24b3e78ba6f, },
        { 0x124c55b2bf346a2e, 0xd103a43d3935e365, 0x057ec4bccb1495fa, 0x00ed10782715baed, },
        { 0xa0acdb33d51a6122, 0x371bc5b61637632c, 0x8e34d701b00d94a1, 0x55a768b27af603fa, },
    }, {
        { 0xa18985d8b26499e1, 0xc1bc62083ceb5996, 0x831cdab35df6f244, 0x649754af630b203c, },
        { 0xa3fbd80eb81e44c9, 0xf7cb98effcfbb70a, 0xf4aa5e2acff94f56, 0x688d5981126d04ae, },
        { 0xb3b29f908cbd0ffe, 0xfa623fb6269a99db, 0xb4b95b97f70e7d34, 0x0d4247a02a21e650, },
    }, {
        { 0x9ea81138312eced9, 0x85e9d4d4fd0e579d, 0x2308ba5023d68e7c, 0x36b2632bf6983f2f, },
        { 0x19a9f3f1f87a7d71, 0x60077d6878b0d72d, 0x8921d76101a59bbf, 0x262e84e2f8407391, },
        { 0xcdc70136d3285253, 0x7fff27a4014124e2, 0xe29ac913db2a2df9, 0x1d3bdc3866338730, },
    }, {
        { 0x583b04bfacad8ea2, 0x29b743e8148be884, 0x2b1e583b0810c5db, 0x2b5449e58eb3bbaa, },
        { 0x5f3a7562eb3dbe47, 0xf7ea38548ebda0b8, 0x00c3e53145747299, 0x1304e9e71627d551, },
        { 0x789814d26adc9cfe, 0x3c1bab3f8b48dd0b, 0xda0fe1fff979c60a, 0x4468de2d7c2dd693, },
    }, {
        { 0x2118278de3bc6748, 0xe71ffd60d0b20ef7, 0xf551be51c67bb198, 0x26a13664d0543d4d, },
        { 0x29522d3b13a339ee, 0x855225506cd89529, 0xdfea3ad4acf4f0f1, 0x49d76bba7942742e, },
        { 0x14fa42338d56e61d, 0x191d3946c351299a, 0x247d576da7adb185, 0x4e1fafe3a8fcedc2, },
    }, {
        { 0x4b413b876b45ea9a, 0xf6c48949fa00ac34, 0x1a1f236281ac11d2, 0x30f3dde7dd980a49, },
        { 0x2d8ba97fbccbc178, 0xaeb46997519fe02d, 0xb84b8fd5b79e5c67, 0x5b91b0fb143a7698, },
        { 0xca900bc53c1c04d5, 0xf169e39fa39b659c, 0x8d2acf0608a5aff2, 0x26416e527c985580, },
    }, {
        { 0x9e6c7295c3366bdf, 0x4808b22db28d0a8f, 0xf7375d06a7082921, 0x072ae4417280cb9d, },
        { 0xad4e95267e73663a, 0xc26e23ff6402bf65, 0x34e1dab6374ecb55, 0x4f1eadf2ac128163, },
        { 0xdf30e1bbc7581fdd, 0xb3e73f2e61a63c23, 0xa750cfd5a16d9c3d, 0x719f9fcf497b7242, },
    }, {
        { 0x15e7053d236a044c, 0x3cddbcb13b8d87e3, 0x519960d2d321a828, 0x4e559a0f0fc5bba4, },
        { 0xfe00e8769c12701c, 0x95dcdc0a039c3b5f, 0xc169454b0c02eb1b, 0x727021d35f87530c, },
        { 0xa571040727df241e, 0xdf45efaab2900d36, 0xfe6edb5c60a69ade, 0x64fcb73007bbc01d, },
    }, {
        { 0x38ef58cc6fd390ca, 0xef786575171a98fc, 0x8850b78fc442d65f, 0x6f34c66d6fd086ef, },
        { 0x93f3cbb43898dc04, 0x0791ffb24307b727, 0xd7bd8096ce34981d, 0x0b598b8e8b849f6d, },
        { 0x11cfc18a0cc2f689, 0x81114607b529ce2a, 0x0a9bc046c00b5940, 0x412128b0b1ac66c8, },
    }, {
        { 0x8b5f16b3a3d84743, 0x809c976320eb9ff7, 0xabc50529df9462c3, 0x1c648ce1c4cc6f52, },
        { 0xf2674254e2eb0307, 0xa835df642da3aa98, 0x377941ce07dd382c, 0x09603b84f60597c6, },
        { 0x1840b2d2c9e55150, 0x0cefa12f390c7082, 0xcb58a4a0b6d2a82f, 0x4ec62e6ee4bf3457, },
    }, {
        { 0xcf391089f571047d, 0xada730b7001d01b2, 0x18871aa4b8daba3d, 0x34d6411bf61ccce8, },
        { 0x0c14a6279f79ef24, 0x1e44af9368276ce3, 0xb72ef4ffaf5e5ceb, 0x72237bad4ac56322, },
        { 0x9e293ed737c0dc1f, 0xe165f9d6befb66db, 0xac1b10fc325dd9f3, 0x3e5ecbe48009d865, },
    }, {
        { 0x59733185e8b81044, 0xe2e029e4108bac60, 0x9d9dedf46908254b, 0x4f8d3adfbfc317ed, },
        { 0x078b733425af0dcc, 0x0472d1d35067ba11, 0x37848b5238191434, 0x543107b9bd371a71, },
        { 0x6b64607cc08d3984, 0xd44f8eb9c89f5b2c, 0x805cd88b77523565, 0x0ce8c056ae0b9ccb, },
    }, {
        { 0xb8b66511b3819ab4, 0x5309b011c418154b, 0x514f7cb560e07caf, 0x47fd085566a3f80b, },
        { 0xba4eadcc76e05caa, 0x9d40e8783a54a171, 0x7c9542ed26da9727, 0x46440c7a2f23cdca, },
        { 0x80102ab0431334cf, 0xcd1b278758c0ce70, 0xda3ccbec0154b9bb, 0x233a3339740e3e26, },
    }, {
        { 0x66b1db8603ff1bb4, 0xeccbfd7b4f427346, 0xf5c0ef08730e1064, 0x6ad4ef7098ac2743, },
        { 0x6b409241a78b10e4, 0xd33139517a431c57, 0xf2423b526d66df43, 0x223ee3935b0b4ec0, },
        { 0xa882277c9067e8cc, 0x4c113d3710a85f5e, 0x77a224c5569ce51f, 0x77b028b623ff2f06, },
    }, {
        { 0x5b209f060d549385, 0xa08976f763c235e0, 0x3d4c47b65a3054fa, 0x60658f6dceec9d5f, },
        { 0x3cc5fe1e06156933, 0x56cb45827c4d221a, 0x913cdcdd6bdf48a7, 0x1ff62032913f6b94, },
        { 0xe5e35afa7095b31b, 0xadebe9a8b528bd84, 0x774ac00f111204f7, 0x73473323af32ba8a, },
    }, {
        { 0x6cc34769f4f67b48, 0xeca62a7c25194e16, 0x1d302a919e71de23, 0x3b64fcf1a0ca9c66, },
        { 0x6c0f0b9c58b1b8c0, 0x3795db5d293bb326, 0x5c8e6a486b411555, 0x29275db665336921, },
        { 0x30c8f7720e8b907b, 0x7dacfaabd18547b0, 0xa786345fffed56e8, 0x7418665a089751c9, },
    }, {
        { 0xcce6bdf4a1dae8e1, 0xd56eb17c9bafbb60, 0x3dd92f761c286590, 0x5e5c02da43f38aed, },
        { 0x0b884a8e366c497d, 0x2d68611c50f8071b, 0x6fede5e94f7d7ad7, 0x1508210fc9190a34, },
        { 0x380903e6a50625cc, 0xf4e3c1232de33381, 0x3ce706816dbd4874, 0x3437de93b66b3a9c, },
    }, {
        { 0xdcb489c7410a634f, 0xb966ca7bc5f43c21, 0xc8151dd3627e497d, 0x7ea1541e444c095b, },
        { 0xb443acb80ee54670, 0x34ad3ff2b8f40c8b, 0xba037921689d8b4c, 0x405952a41a082a71, },
        { 0x56b8d8895142abc6, 0x86ed618f540de993, 0x72a0cd17cc7fbf1e, 0x044745ca468b4e24, },
    }, {
        { 0x58fa4ebef5bb9f3b, 0x1977b7ddb99c8573, 0x07ebedbb81d0e028, 0x6491ba77c9609c12, },
        { 0x2fbd585aff8604b6, 0x9fe55125c0b44e60, 0x9ba69f90a3aeac8a, 0x05393a0615be487d, },
        { 0x67cd556fbb6b377e, 0xa89d57ff17735944, 0x38db21670c137883, 0x3b77872f90bc96cb, },
    }, {
        { 0x96eebffb305b2f51, 0xd3f938ad889596b8, 0xf0f52dc746d5dd25, 0x57968290bb3a0095, },
        { 0x4637974e8c58aedc, 0xb9ef22fbabf041a4, 0xe185d956e980718a, 0x2f1b78fab143a8a6, },
        { 0xf71ab8430a20e101, 0xf393658d24f0ec47, 0xcf7509a86ee2eed1, 0x7dc43e35dc2aa3e1, },
    }, {
        { 0x39b3240a8e260e07, 0x852074f462871ac8, 0xe2b30d7ba44e0939, 0x46024d80a0abaae8, },
        { 0xb2b999b9d3c313bf, 0x2d5abd526e51ebb2, 0x90bbdbbbee0e550f, 0x68da131c98ccce19, },
        { 0xec23de8e54a41b0a, 0xcf7b9637b782780d, 0x3cfd7049939e311e, 0x597a78170ea55972, },
    }, {
        { 0xbc045fde06baafda, 0x47739ae0fb5a4416, 0x7c8463d85d2ee4a6, 0x2203eac504374ad2, },
        { 0x0eff92fea05c01e7, 0x75f7039c895fb7cc, 0x8d6a82179404e0f8, 0x573560a58a9a737a, },
        { 0x86782d1959dd3897, 0x6012366770d83750, 0x9b25d884fce034c3, 0x3e1460c501e23583, },
    }, {
        { 0x0532f811016c8b1c, 0x2e95ed3fd136ba45, 0x63a63935e61135ab, 0x4528b81a998ddb46, },
        { 0xbb03507cb2e76cb2, 0x4974972d6fa0507e, 0xac4ef385909195f3, 0x06057d24584f70a3, },
        { 0x54a5662f8c90f978, 0xcbc6ccf52f7b10c4, 0x448b019b865db58b, 0x42420ee395f99ded, },
    }, {
        { 0x30baf789b5714ceb, 0xb8aa475f974fd937, 0x433fa7650f057b03, 0x1f64562f4e4b8810, },
        { 0x563ee27295c2fd9c, 0x5dd036c4a9b8dc68, 0x7c530ef9036007ba, 0x39f11784dbd2061b, },
        { 0xab6c450037138fab, 0xc27f331534f6dadc, 0x0210194c8f069781, 0x50e316319e2f307a, },
    }, {
        { 0x472ea0bcc3b0007b, 0xe17ca69e45eb13c5, 0x6278a077819bfa26, 0x09ffc1dab21b3d0a, },
        { 0x335f8d6c8db01743, 0x2b131c1c476c43b0, 0x6fc7496292ddc403, 0x4fb1e3e11b87720c, },
        { 0xe7355cbc93b99dd4, 0xd1c3753568edcaf2, 0xf8c2a48a113c1f1a, 0x78a7a0577d5eb4b0, },
    }, {
        { 0xbe51146960456cbb, 0x1f30b9fd8d5d8c2e, 0xb87b3385b4c04589, 0x554247d44351534c, },
        { 0x101625e255283c9e, 0xfad2395dc1e61fb8, 0xc35340596b7d5564, 0x521864bf96e1d1f2, },
        { 0xf4ffc63ab9332cc5, 0x918219b2a52bd4a8, 0xb9da028e841911c2, 0x6b1fc8bcc36c5963, },
    }, {
        { 0xc69c8cbaff5efdcf, 0x2911c7018ef90996, 0xc9245fa686929cd5, 0x680dfc97403c2507, },
        { 0x86c5ee26c62e53b7, 0x5f38e8d3217cc5d7, 0xbe8f603d3be42bce, 0x1fa5774d0c2b9113, },
        { 0x89ecbffb2670beed, 0x02950e5d6138312f, 0xd496161120860ea8, 0x1f301f21e64c338f, },
    }, {
        { 0x39d004c975e3fb78, 0x18b6e86529b114e0, 0x38ebe25fd16c3194, 0x607d23449e79f5cf, },
        { 0xce240f2c2338841d, 0xf0eccd758631e5d1, 0xa57ee5f1f18a552d, 0x7ab87bc81c9b3233, },
        { 0x7e7f71501b17760e, 0xb6632a86e158a26e, 0x39ea5398be0ac1d5, 0x68a56d9d503ba911, },
    }, {
        { 0x23cb640ce76fdb4d, 0x47e311133f0908fe, 0xcdcbac03d0b34bd7, 0x7f82dc07837c25b8, },
        { 0xd8e96e5271a8c4db, 0x23294a507cfdb8c9, 0x469606658af236c6, 0x399a2e1dba025941, },
        { 0xf4286f373a40d717, 0x1b14d057a6c7ab20, 0x83a1c82252f67cca, 0x7cc9285e1b1ca416, },
    }, {
        { 0xffe05716c8a5bc9b, 0x066a99e4465e9f43, 0x5f1f7c8b4a91d729, 0x24102bd91c435b3b, },
        { 0xf20b8f9ef4ca0f68, 0xdc05f306e0aef35b, 0xb4666f042d531e9a, 0x383ff5d31b535ee9, },
        { 0x84dcd7cb442b19aa, 0x4118db5fe6b5984d, 0x4fe5dd4621a7dec1, 0x1cc49092d86ef74f, },
    }, {
        { 0x8b29913ab44c3a2f, 0xcd694c5daa52080f, 0x18f938895add88d4, 0x63eaaac5901bc6a7, },
        { 0x1d1d5492b5715b8c, 0xa4200d1dd839a828, 0xbcac48242b183d76, 0x0cbb60c54c296927, },
        { 0x80fbcdfe4bbeb387, 0x2a83408cec5778e7, 0x1f78d31b34339af7, 0x3f9d64c8c8a8ef0d, },
    }, {
        { 0x9dceebbe811a5a5a, 0xa326670692764fee, 0x90359a3c3caf98aa, 0x2b83d213ff96a438, },
        { 0xfacdc91fd47173ab, 0xadea7b2c5d2f974b, 0xc13603d06fdf8e35, 0x30435392c639ef42, },
        { 0x1263d1ac8cc47e2c, 0xe21c9e3e0ba468d7, 0x537a318ff0204559, 0x10aa737e2fb95ac5, },
    }, {
        { 0x85b007aa2006de8f, 0xbde46c9df6bb993c, 0x7bc28481feb58e38, 0x78b9fdc36212a44a, },
        { 0xa5ea6c05decb7244, 0xd24015a74e0d68b8, 0x546c84b69c10e916, 0x2b5447c67bac682c, },
        { 0x7365f03271687f56, 0x8289c896d4e60ea8, 0x419a7cb6e6a09396, 0x3dd7f389d6400ef0, },
    }, {
        { 0xcc2fd8e0a70e980d, 0x6ddf70298381b94f, 0xf50664b1480c7f27, 0x49058f81231ac201, },
        { 0x068acf532b002594, 0x2db03c1651bde407, 0xae564bdba06cb4bd, 0x4cef9672345e692a, },
        { 0x83f6accf64ef121e, 0xb662c5cd03c67508, 0x0096e0895e8ec922, 0x651f80626365f801, },
    }, {
        { 0x16527736446a6929, 0x065115b56a71196e, 0xdb7ce2f7f29c4eab, 0x21e2ae1677cce053, },
        { 0x07f1b268163b9183, 0x53a981088e5270c2, 0x6c820af9fbaa6979, 0x5962bc13fa55e6e9, },
        { 0x287e4e5ce021520d, 0x6863c66494c883cd, 0x58163f9d5790868f, 0x753907ec7e31c55a, },
    }, {
        { 0x9aa6b27c7b784751, 0x324bec2c478ebaa4, 0xe9ed08e66aec068f, 0x2201f3381289ea7f, },
        { 0x88b6f93f1d82d9e5, 0xc753db5538af97b8, 0x1bf3b1612390c879, 0x138e33b33633ed25, },
        { 0xb6b55f377b9cd213, 0x2898d18f20636109, 0x0b5507eec10d1935, 0x09f8a06ae7915b1c, },
    }, {
        { 0x2d49686ea156bee1, 0x2324e89c1c1da03b, 0xd04411d1ac203d50, 0x4c6766231ba09051, },
        { 0xd3afc637aad02187, 0xdbd9df224a427480, 0xa16675060a76a644, 0x4d4f7873219197e9, },
        { 0x62ccf2c67c9fe358, 0x9247d76e62498039, 0xf38b6cd8d7566abe, 0x75b11cf8c416d3ce, },
    }, {
        { 0x263738aed35ffb34, 0x832a1be508f9d4ce, 0x85a8a919b79ab05b, 0x249a81e3a2783d74, },
        { 0x715bb0531829704f, 0x1ec6ac53543da8b4, 0xfda431a312478e5e, 0x499021e4496d50a6, },
        { 0x899a7ee9548cdcb6, 0x6c1dd7b622e94b30, 0xc83ac3b4cc3bc5ca, 0x7058121ddd006024, },
    }, {
        { 0x655ade93d0a8e8da, 0x54028bc0f3b0b528, 0xfeeb9e377068729d, 0x2f6d035763d404ab, },
        { 0x319e7e4057911cd7, 0x48cd05564a053a89, 0x8db0f64a1711b185, 0x23f42f873cc73fd4, },
        { 0x0e966846fd5fee8e, 0x5e0bb86234d2c0ae, 0x60d4f63b80fd5963, 0x68dd930ad8e4dbbc, },
    }, {
        { 0x9fc17ecb49349375, 0x5e732ffb4158c55c, 0x8759bfabfbdc6ded, 0x0a052f827b86639d, },
        { 0xe04cc8d1777962b3, 0x8835c29563464020, 0x7d2e941f10f7a751, 0x50400cc16ee28e97, },
        { 0x654a5e35e09b6bb2, 0x63a31e9c7997775e, 0x78190eb83416c504, 0x304dbdf52bb6f3a0, },
    }, {
        { 0xc0f200e99fca639c, 0x54d06b3235f7f0e3, 0xed9121bc363b7d83, 0x3e5ea2d94acdabf0, },
        { 0x3f317c5595849740, 0x639a7ed4bd56a2f2, 0x1ccaa5f1529c3a7a, 0x5652c4f50989b5b0, },
        { 0x1bb512bbb99e3a9f, 0x4efc5e35bec8021c, 0x1beede8c687c6e35, 0x6d7f7cd4fc28959c, },
    }, {
        { 0x37f602fed6ae865a, 0x308fbca25644c995, 0x0ecc084c77eaeae2, 0x7db8fe19b4cdc74f, },
        { 0x2641e3591377add9, 0x3f7e961cb1a41953, 0x968afee2401f2f9f, 0x2e92a64e5a37bbef, },
        { 0xf20f479793acc7fc, 0xc548575540bab1e7, 0x7c83e40e7f03710d, 0x69f8fd426f44a6ee, },
    }, {
        { 0xf381ad3b56a50dda, 0xd6de8923f2e4a6fa, 0xb375936264d7dcdc, 0x3ead5cd67c2e641e, },
        { 0x910afc9f845b8d5c, 0x714ac6e6a1a43909, 0x80f27d9306b29a0b, 0x6532988a1e41f09a, },
        { 0x552c67dfc88ae13c, 0x8f591c65bf55730f, 0x5b519e2c2fe7a905, 0x6e7f879a77e7a32f, },
    }, {
        { 0x70e895756d2128f5, 0x11b2e65d113e2943, 0xcc14522f502a27e3, 0x0e6258c55223d989, },
        { 0x5f93903582ec85ae, 0xad15cf958535fd78, 0x756cdb90f786332f, 0x2d18d5ed5da65e0c, },
        { 0x8ad6b63e9bdfc49a, 0xa1b3286ae7f7c76d, 0xb0e423d856d495f6, 0x1461d531acd6fbc6, },
    }, {
        { 0x277b0ff7136b2904, 0x6f94fdaea8d39066, 0xe9a20fdc56101bd9, 0x4226ab45b892705e, },
        { 0x9770dcb15f57047b, 0xaa2f99a9c5cf4123, 0xca99da38bbb7f32c, 0x5b74a5196a654557, },
        { 0x2c3aa0d1494db970, 0xf8fdb612af7ca799, 0xd6e75deec7b28ae3, 0x0ca5047c982b6b0b, },
    }, {
        { 0xa096097e5d5ab1f3, 0x55e6ac6b30d14c4d, 0x0e787b5b40b47347, 0x27fdeb5f7172e37e, },
        { 0x5f328ad303373c3f, 0x2e79ad7e3b3038a6, 0x39dc8f94de222f9d, 0x34f4b6d11b18555a, },
        { 0xfbd02a20fca72fba, 0x2827acb24fcadb53, 0xebc9a4e0ee77588c, 0x7afcf628ed302543, },
    }, {
        { 0xb85a436c5e1a397a, 0xe51b3d09abfd5dc6, 0x09509cb9351e57a6, 0x25a1c873b3c8f3db, },
        { 0x86f2bebe23fd941e, 0xf568f06cfc6e6286, 0x736ab7bc2e4d97a3, 0x682c7169999bb98d, },
        { 0x0b156f550d5c4162, 0x2482776f2754d6c2, 0xc2d1f94dbdc504ad, 0x1d018e5ed5c203ca, },
    }, {
        { 0x43bf773529ed5b5a, 0x5edd0e519be0d0da, 0x6d4264296da65344, 0x34a7969bc489c0ed, },
        { 0x1643ab04cd5d792c, 0xd91d4d117627e597, 0x42679317c4545cea, 0x4d661e38268698cc, },
        { 0xe6ffba43d7ca94c6, 0xfd395c385bffbce7, 0x2790702826310616, 0x5e78150a96f73c17, },
    }, {
        { 0x763966fcfa433f7c, 0xb90424c06903b07e, 0x2a9f9b88ce43e341, 0x55d4c871a8eda21a, },
        { 0x4860d96d88047ffe, 0xd5a0d55097ccfdab, 0x05db5a13768e3fda, 0x3f33f0265ff6230e, },
        { 0x48698d99183aa8b6, 0x8beba46ae3871279, 0x40f762bc7767ee98, 0x26bd664e6913b8db, },
    }, {
        { 0x6e61c09ccbf721f9, 0xaf380d23d47ef6e1, 0x82784ae2338eab76, 0x31816ac9167f4fde, },
        { 0xdee539f6c6bdb362, 0x591d41449960c22f, 0xa5b1e10da3a085c7, 0x2ef298a42fdfa6fe, },
        { 0xea4492c9a6b817c5, 0x72cd8c86961b2033, 0x83faba6088f3a779, 0x66012e6abbb5e2c7, },
    }, {
        { 0x3289bb632f14c065, 0xd4a81cfe95e2f595, 0xb616cbdeabf8b5b3, 0x5c907f1f7a49feff, },
        { 0x17b5d1171a7c3ade, 0xede226728f5fe973, 0x184c3adf54a3c3e3, 0x40ddedbb642dfdca, },
        { 0x5f43c05179bc659f, 0x5324fae02ebe945f, 0xc35fbcbe6798eb6a, 0x065e8bc22f0245e7, },
    }, {
        { 0xa66dcc9dc80c1ac0, 0x97a05cf41b38a436, 0xa7ebf3be95dbd7c6, 0x7da0b8f68d7e7dab, },
        { 0xef782014385675a6, 0xa2649f30aafda9e8, 0x4cd1eb505cdfa8cb, 0x46115aba1d4dc0b3, },
        { 0xd40f1953c3b5da76, 0x1dac6f7321119e9b, 0x03cc6021feb25960, 0x5a5f887e83674b4b, },
    }, {
        { 0x0a8d60180ca2c1f4, 0x815eb0dbcc68df40, 0xd7e67a47b82f4e99, 0x45a02890607f15c0, },
        { 0xfef366d1fd41f184, 0x8b694a1101cfe11e, 0x4b39e15e0150a74d, 0x4013f03d6ad351ba, },
        { 0xbd0282dc6ee065cc, 0x36b994fd224ae646, 0x534e9ad8febce874, 0x482255c1d9f06e4f, },
    }, {
        { 0x1e6b0a5605928a7a, 0x1d5b915125db60c2, 0xfa8792c0cfc4755b, 0x2086f7fe55427d2e, },
        { 0x0204ac4e8d2ca3ea, 0x299bd9aaea1927c7, 0x0e6883eaf325cdc0, 0x62e9884065b65b44, },
        { 0x7775a78977faabdb, 0x054b62586bbb5798, 0x210a7ebc9568ba16, 0x1dad6e30e60ab511, },
    }, {
        { 0xef459f08628ed648, 0x57a3c57b90e1486a, 0x54eb20d2eb699406, 0x298e3ae27c2a796f, },
        { 0xe221b34a080f1721, 0x394438a7c1128e4f, 0xc370ce8b7d18ecfe, 0x55c5e6a7bf53150b, },
        { 0xfefb6e91b03c4ab0, 0x54b38d6f20fee7c4, 0x395a19af3c3d04ac, 0x79ddf3a97a8b9170, },
    }, {
        { 0x3c03eacf71cef800, 0x90367544ca8afebb, 0x383fea286a29c477, 0x4e8593b0bc655462, },
        { 0x12de114aa3e5638c, 0xba2a4aa929c4f20d, 0x56b0d29d7b8b13a3, 0x6bb91a497b9b7944, },
        { 0x2a49e646c5e7d206, 0xb13ef9cd9263c445, 0x50ab6ce8edab529e, 0x20cf7d79b0ebe39b, },
    }, {
        { 0xcbd28f4e8ae75c48, 0x3cde029144000b60, 0x373bb9c898bc2170, 0x7c1188539f570886, },
        { 0x7db4939df0fe7dca, 0xf50eb90fcba951ce, 0x098be61c357e1d1d, 0x023562378899469d, },
        { 0x20f6effae15a4c03, 0x2f470a943c778e05, 0x79f50a03fc99de67, 0x38d20188d1061483, },
    }, {
        { 0x05b6b816118ca612, 0x762934232acb38ee, 0x3addab78d4c85d20, 0x2872e796f26be705, },
        { 0x0f2d95b40ce21195, 0x14029e27e160884e, 0x3dcd49f76a05bae9, 0x23b7de357beb7ca6, },
        { 0x78d00cbac0c99a90, 0xc48a93f4964dd704, 0xa010301c61b7da2a, 0x4e307a033f0da51e, },
    }, {
        { 0xe39068107d90930b, 0xc4d5c0190926081a, 0xeac881987eee2811, 0x369e7a037166a974, },
        { 0x4eb04d2ef32a7aa7, 0xc333e92623486432, 0x67f2859d5bb4bbd2, 0x340611ce362db995, },
        { 0xdde8f5ca69c77c9c, 0x9b8ae19d5612a945, 0xf989fefb3b0b7046, 0x3dd66b1a2fa594d5, },
    }, {
        { 0x4daa10a21f67006f, 0x1bbd786ed2cac6bb, 0x5668e51b34ef1aaa, 0x76ce237084abc92e, },
        { 0x775b57c56942f699, 0x44dd2b3c62952644, 0xccc2e090d6b3ad0d, 0x3448386f929d4430, },
        { 0x2f971db350b99804, 0x1cee220f59461c34, 0xf0ae8f967b5586d2, 0x742ca51fcbb727f9, },
    }, {
        { 0xa97d2d2740a109e2, 0xef72a5640c5ade76, 0xb3153e4d47f89224, 0x0757179ea9c2453c, },
        { 0x0d4300bd241a5382, 0xa0309e325e919492, 0x8d4b6e261c896343, 0x0eae2d6ed576b0c2, },
        { 0xd6dfba4ecb862b52, 0x96868e971106acee, 0x2564160686e7a06e, 0x6257151a9fae1824, },
    }, {
        { 0x6c99d2ac2182d28e, 0x047a4911efaa861b, 0x10593c99971aee9d, 0x1720282987cd88eb, },
        { 0xd5393fb2f6293845, 0xe297eef3a3cbefd0, 0x5054c9530e6d7860, 0x4849b391f5dac84d, },
        { 0xcc291020cb38dba7, 0xc3af4d407e4405f2, 0xe74b6343e9008947, 0x0ae477139c6676fe, },
    }, {
        { 0xaddc19fb808fc0c7, 0x833f2820b97f1833, 0x3cae0e4ed3527f9a, 0x38ae17a2193e98cf, },
        { 0x78352d6f66c546d6, 0x6acef560f767b9b0, 0x8bfcb4a47943d1b3, 0x58f6a1757ce7a57f, },
        { 0x9dcc21a4418194a5, 0xdc49b44e460e3a36, 0x7091ca3d197dd7e2, 0x3bb4dfbad8e367da, },
    }, {
        { 0x12c0ff19b44209af, 0x2bc53cce0b46dc9c, 0x4a6a7c7119ca4bcc, 0x5c11757f98f3d04a, },
        { 0x13ed021320cc12b6, 0x3adad589c82aef76, 0x8dbac08502ad8434, 0x1314f390e5d56909, },
        { 0x26145f55ebd6c1d7, 0xa18896b47fe3f892, 0xe8de0d1fdd90ab3a, 0x3d6c810d2453bd12, },
    }, {
        { 0x7d53625f21980c7f, 0xeb61d5c748e60711, 0x733a560c6cf4da80, 0x1cca65ea1102454d, },
        { 0x0e0c59599a758d78, 0x1c89bbb38649a297, 0xb1ab7d29c29ac6e3, 0x326436e80475fa29, },
        { 0xd0f89edb123b09e6, 0xcf0fd49d9fe8d47b, 0x4a9d68678c24f02b, 0x3968df7296df0a45, },
    }, {
        { 0xe3bc39ca2bca971f, 0x19e76df79bdcf7af, 0xb37414952a65a57c, 0x7e6fdf274e5a0fa4, },
        { 0x1d8e23b0d9f7828f, 0x0da52d29891fdd75, 0x3895641b20971626, 0x4c5b05a206756e8d, },
        { 0x3d16f1309a392a60, 0xd588033f814139d4, 0xaa29b40c65382c75, 0x10657167e844e19f, },
    }, {
        { 0x45df187593b015f5, 0xa4b4133641c1caaf, 0x7525f9cda362d173, 0x096b4489ce1ee76c, },
        { 0x466c16c7f22f1fda, 0xded7df054c8a69bd, 0xaaa4793350d4be1b, 0x5419a8861aa0a635, },
        { 0x4922ca3f1f50992c, 0xecfd2a631c926d3e, 0x2a958b98783afb74, 0x5245b5e710a0a934, },
    }, {
        { 0x23e811ad0e6315df, 0x0b650d05e2aeb290, 0xb7ba0f59a75d586c, 0x043eedd45e1f4dee, },
        { 0xf6c147f2c7073217, 0xc651b919f3afd20c, 0x258fdbfd7041f802, 0x173c4fa94f45073e, },
        { 0x3d71ea60928df9c4, 0x5b7e78063373562d, 0xd9b0514ca29552b2, 0x1e2a7024993cc472, },
    }, {
        { 0x601a0fbcd45c811f, 0x24b7bc7d92ec0803, 0xa0cae62b17d2407f, 0x5fcb43ee06225b26, },
        { 0x310509b93509fba4, 0x0d8db37605631b75, 0x97deccba52401c87, 0x044649f411b2e773, },
        { 0x0c0d24ad9598215f, 0x1b7f9026cc36628c, 0x338e2f557016dcea, 0x0c8a1bfa5cc0e58f, },
    }, {
        { 0xad0e8d3e5a498632, 0xbc5472d16664c2fd, 0xd88564919b093759, 0x0ec4dce31c6a902e, },
        { 0x0a7fb14ca53b3119, 0xafcfd591ae58286f, 0x75c643f5380edd5f, 0x76023803c4bb9c9b, },
        { 0x902b8965698035b7, 0x7948c088d5d7e42a, 0xde7471dfcb0ee751, 0x0976a7d9cc9cb0b9, },
    }, {
        { 0xed8930678c69bf9a, 0xf6277cb9395a26fa, 0xb1b45ae5b2423bc3, 0x38d6fea173b49b63, },
        { 0x5bb02be436c79869, 0x0fe7b60213405daa, 0x78e9e6add83397cf, 0x0e146c2835ce3966, },
        { 0x91a0862a4b0725f5, 0x1e29935460233a46, 0x770c6408331cd827, 0x0c3b6e0b04fa705a, },
    }, {
        { 0x8de703b5681d104c, 0x3d2f7a591263cb45, 0xae710c171ce56c63, 0x6b857c7efcc3e6ca, },
        { 0x79d256b48b2801c0, 0x7e9fbeac3c400fc4, 0xa751ab1d4733ba41, 0x09de2bf5dd418aca, },
        { 0x3bf10ff3eff0687f, 0x5ebaea34f1e37ba2, 0xe49e61261d66034d, 0x5b466e2ac3b242ca, },
    }, {
        { 0x137eeb6747fbb842, 0x79df5c7560811a8b, 0x5a2ba76f71f8c89a, 0x09952a563bc8ffc2, },
        { 0xa2a8cb4bdc7ef83c, 0x96b5c6fa5f93c226, 0xd4ebeb1b0664e3a5, 0x409b4adce5c6cf2f, },
        { 0x44d53db9834350c4, 0x89299305a5f505b4, 0xfb22faa25949ff2f, 0x69b968a704657d64, },
    }, {
        { 0x6c9ce663f3b83f1f, 0xf12e893f9ce2ae57, 0x04efc97232c4c3d0, 0x5f58fb4de69732b9, },
        { 0x04676e045348092f, 0x0eacbd686fb6c492, 0x56091566a6d70147, 0x62d7bfa65aca06df, },
        { 0x3963d9e7e2407206, 0x0eead76f785795e1, 0xe134918ea63bf4fb, 0x65e052afd3447ff9, },
    }, {
        { 0x56626c089a105205, 0x6a68847e4ca91191, 0xacc899c90aac0cd7, 0x38ffe58e2fe2a51f, },
        { 0x2da208317d397715, 0xb3a08734716771bb, 0x5e5eaf8cdab3b832, 0x735fdc7202f0feda, },
        { 0x95242682843bcfe6, 0x6cd54f0a3c2fc05f, 0x10835128a8ea74aa, 0x019812eb2ee8a03c, },
    }, {
        { 0x77257b3a43495738, 0xb8ea4be4c51d872e, 0x415abaf7703993bf, 0x5eb64903987df031, },
        { 0x5780874c187c9965, 0xa1c3f5f59132ca0c, 0x5ae601625e6afd30, 0x74ebfe7839c2a855, },
        { 0x8a1c51444a608a8b, 0xec97a2f38857dbb1, 0xa4ea18a3793f5123, 0x42e2a9efe1f42817, },
    }, {
        { 0xb5ce6f2cd77e277e, 0x3277a6da72cf0ca0, 0xa39dd9ba4815441d, 0x2cae3babe44f7d25, },
        { 0x5b6e0a466eba994d, 0x1c1e716839936244, 0x648b80ae9456f39e, 0x37298a1245e72dc7, },
        { 0x2df48b8b1b1a6a67, 0x66613c4e625b13bb, 0x613d15308cfeb944, 0x78bf5dc2599e09d5, },
    }, {
        { 0x508d09740104b715, 0x40da1307f497cb8a, 0x4dbff163ed9c492b, 0x4bd7af7a23943739, },
        { 0x28d96765409d02b6, 0x2666f83f945a433a, 0x4b644b7fa708399d, 0x1158fbb4566fb51d, },
        { 0x4350a962b0ccd019, 0xc3308270a9863f9b, 0xd103fbaae057038a, 0x745e1908b3015d8b, },
    }, {
        { 0x1d965a13ae62fe0a, 0xe4237e736ae6b9c0, 0xd41a5b18a13a16eb, 0x3421a6864e6c256e, },
        { 0xe07ca9f743ddca02, 0x0c61c72673b5945d, 0xcbc8273103365ca0, 0x33b2aa2f2e1b8b3f, },
        { 0xa4cac55070fce4d1, 0x214094084b9b5f95, 0x1f5c791e2a1d42c0, 0x1e7958cf4764bf69, },
    }, {
        { 0xba9f9510c787a659, 0xd788b54303ed5493, 0x683b5fc0b0a5569b, 0x797c8e25c0fe3929, },
        { 0xd80de6381330964e, 0x420fd04ae21238be, 0x07a1b3e3a052a0b1, 0x61e8f308b9f68790, },
        { 0x31bd7f2bf1396f65, 0xb111846f2419860b, 0x87b636d141030b36, 0x30e41fceea50d463, },
    }, {
        { 0x49d07989119b232b, 0x7af92508352b78b8, 0x2998127625d656c0, 0x23b4e084c20dbf1f, },
        { 0x3f963735eef5658f, 0x9ee0debbcaf6cd78, 0xfead3bc0487f1590, 0x70330bfea2da3ba8, },
        { 0x68618aeebb29c84b, 0x886c8d3d44f6b1ae, 0x55d8fbb738af923f, 0x45838d3b0cd689de, },
    }, {
        { 0x84ac4fe361295fb1, 0x1cdfd6757bed5a8c, 0xe99428bca88212db, 0x373e0b2a34eb0afa, },
        { 0x6ae452dfcbbd642c, 0xe5426167e72f33e7, 0x4434cd39c0c30398, 0x33004dad16ac6121, },
        { 0x5b436c32707bdeb2, 0x73f4492a54093925, 0x20f34a2534121367, 0x3c15a2ba625da7c5, },
    }, {
        { 0x86a977f58a61c57e, 0x0e4069d5bf8613f0, 0xe24a7e40b0d80789, 0x0b1cec760b4f5bf0, },
        { 0xa413af861dc892a0, 0x4cf5cb841bfcb0d1, 0x25c0f20ec0dc4fdc, 0x0d65fca550b64c66, },
        { 0xeb413f935b05e655, 0xe36b4c6f3ef0cc03, 0xb256d88e0df058cd, 0x07adcba113722755, },
    }, {
        { 0x79d7af892f4951df, 0x88a06eb5afa6eb6d, 0x5831563c24c568e4, 0x66707a19aeb26ced, },
        { 0x6c1310e2747396e1, 0x0ee400624cede7ce, 0x9a3d76022d45e29b, 0x204cbd043f0cf648, },
        { 0x2b2454e8379e42ae, 0x90c8758cf12c7764, 0x6c5fd2cf65e7c6e6, 0x289e7f37cf80a50e, },
    }, {
        { 0xa16c84264708cca4, 0xb6980348e436aa70, 0x4f2683d04600a6e9, 0x47095bed3e6f41eb, },
        { 0x9a971f3dea0ec072, 0x787a1c8acff14f17, 0xe1f6197c49aa1aff, 0x24705841bc0551f0, },
        { 0x659aeeef135be2b4, 0x6ccdad89fb75b6dc, 0x014f9d4f6b16e15c, 0x1e7b503e03e29b5c, },
    }, {
        { 0xbd217b563fd31b95, 0xfd09fb024342853b, 0x52ecc0d04901931a, 0x531f7220b46fbad3, },
        { 0xb0e55e3cc00ba975, 0x8937c6c0c8225fb0, 0x00846f52fd60e6fe, 0x159452d6bf5ab08c, },
        { 0xc98de4465688f03f, 0x89ff073b0e6107ea, 0x3648578e9bf173ac, 0x3789df6878073a29, },
    }, {
        { 0xb11d5529e03e12f7, 0x75eb970fddf46f16, 0xbe1f0518f885f334, 0x2d27012b68e343c7, },
        { 0xe90965d7a91ecffd, 0x63bbdd908234d517, 0x9380df5651719d88, 0x2d096724c8f44a15, },
        { 0x7267e789e2cdc6d4, 0x714f8be23672f87a, 0x55cd7495ea9b45d5, 0x1f9f97449943edee, },
    }, {
        { 0x5fcb17c89b58521c, 0x7dcda9176645c49c, 0x3136efeb5917bd77, 0x7d5ef8a93673c3c1, },
        { 0x7f6d65018666dd69, 0xb013d32e2780c2ca, 0xa607d0b13789882b, 0x6e693bd4b5554a65, },
        { 0x6337649c513dd959, 0xf3489531a7717ea3, 0x0bbb78c72e946717, 0x3357ca52814bcce5, },
    }, {
        { 0x7e720b523ff58085, 0x8fb70c10ec6a424a, 0x925fa8df3fcb64a8, 0x437e113dacc6ce63, },
        { 0x88c0114b01fa6a5d, 0xc168d1cb19004a46, 0x4d148f642e00cbc1, 0x58c5d1d92441151b, },
        { 0xd5a8d6624f3f6b4c, 0xbcebdd2d951aab88, 0xadfd5e27c2603ace, 0x70ea32335dd3e15f, },
    }, {
        { 0x9514265da9684d92, 0x55419cf8e74d2bbd, 0x34451d1d2948c133, 0x2420b420eb44b9ad, },
        { 0x7b9f2ec17d7d252d, 0x49b2eb299a152ee8, 0x11a954c54e9359af, 0x2f4d538eb369bb94, },
        { 0x2ca37418e2bb55bb, 0x5819d1000bf6b294, 0x9e192e132832f057, 0x3b99ef9d2281b8a2, },
    }, {
        { 0x792e12e456f4b677, 0xe0e7dd510bd62435, 0x53a50096621ecaac, 0x651b2093bf8dd44f, },
        { 0x7c34ea995bb90b64, 0xed6b5c4f498827dd, 0xbd2c09fb7e0cac91, 0x2e87c003fc3f6c76, },
        { 0xb9a05483cb39f86d, 0x396caf957e208eaf, 0xa49f69b492d84129, 0x5bd3cfdcf0aed850, },
    }, {
        { 0xa81c60891106f4c8, 0xcea547c648eac4c5, 0xbabc5e19b65a4461, 0x1edd6ba7914405bc, },
        { 0x512e6fb2b7a10a14, 0xe5a5ed1991fde918, 0x9c7ae07230352098, 0x7480322fa86bcbaa, },
        { 0x67d298b8b0b208c5, 0x16c79b83fc05876a, 0xf2d11ce71b09584b, 0x718679956d0f2822, },
    }, {
        { 0x5bf82d6078a0ed66, 0x897f84ec9b62e796, 0xc93d0b6f351ebddb, 0x2a530e07978dc6a1, },
        { 0x25df69afcc5b4c55, 0xd299fdac4b1e7b43, 0xbcb7e7f4dc8026eb, 0x0f51fee0f671a838, },
        { 0xe3abc271e051a08a, 0x5e9c303c43549c50, 0x971ee5c96e751049, 0x6af8ce480990fdde, },
    }, {
        { 0x0d2ec63e85072941, 0xef34d76f7629e4dc, 0x4d299fad3eeba243, 0x1c0c8ba81d0e272f, },
        { 0x49ac8b607de8947f, 0xa598dabfe6b60d1a, 0x7fcd7397262bf2c7, 0x014176708b05c98e, },
        { 0x390985788786e66e, 0xb2abcbfda155b485, 0x1abbee0e9c793f7e, 0x493edfaa969224b2, },
    }, {
        { 0x249e0a047ad767f2, 0x1311c2f6a3912f40, 0x897edf4f3f0934de, 0x0109d1cc1f5f9664, },
        { 0xc6db180f36fd7583, 0xd9504d7c6a3f0d70, 0xe317b1c26b2dc20e, 0x0331dc461e0c2804, },
        { 0x213c7da0cf2bab83, 0x72c760fd56f00288, 0x4e13d9f2044ff105, 0x54c1c1f159ea9372, },
    }, {
        { 0xd9439624ccff195a, 0xa69a217b04aa2dc2, 0x7b58138177cdebaa, 0x727062c36c36f14d, },
        { 0xe38363b2f156d356, 0xe2fca0ec1a7a6c40, 0xf3236173336d9cdd, 0x50f8528554708983, },
        { 0xf7067e62f65be459, 0xe22ab98afbb87c7c, 0xcd2abf428127a332, 0x474f259f41243874, },
    }, {
        { 0xbc2d1e80b082f07e, 0x368cea6ea25756bf, 0x53a67cb57af12c86, 0x17186ee8b385ef92, },
        { 0x3a022d647895403b, 0xdf7178fb3eafd6d5, 0x1bcfc70ac50d54fc, 0x36e9b22d4569f3b1, },
        { 0xbb82cdb2ee005d40, 0xb99d5b3901db4779, 0x291f080baef8c235, 0x5326c57b4e3b4d18, },
    }, {
        { 0x0bfa2499f01d5b98, 0xddb2533f8ec91d4e, 0x7168d7e71bbf5701, 0x7d9e5af9c9fafe4a, },
        { 0x7ab3ab447310593b, 0xafc3fa66ae863030, 0x45e5ddfbdbd73272, 0x65401219c668ccec, },
        { 0xa57c099f9b773f39, 0xf8f0223b7940844f, 0xea7c8a4683aa97c0, 0x6860dde450335d21, },
    }, {
        { 0x3f2358cd8a54d845, 0x5e1cbb1b2203028a, 0x63f1d82211fd3994, 0x37401e654308d8ab, },
        { 0x1e9adce407b41d27, 0x343ae381880a16da, 0x49b82be44fc4494a, 0x6fd66d342726aec1, },
        { 0xeca260107cfe6583, 0xf981d58463921d49, 0x1b69d8d3bb516bf2, 0x461c28c6a669e323, },
    }, {
        { 0xcd69db0a01c923e9, 0xc095fccbb2ec7a51, 0x6bd7c77815518043, 0x0da33c93b6bde853, },
        { 0xf3bffb0a349903cd, 0x0362f6bd664c68f5, 0xb75c72e17d7ace8e, 0x12813206c26de303, },
        { 0x8e82eb864d21fb37, 0x784e6aa7d706fca0, 0xce5c50197986f343, 0x749a29cb86222c9a, },
    }, {
        { 0x5698f7bef14ee1f4, 0xaf410069a93c6b27, 0xc3b2308b2f72e02d, 0x5853c9f30697a906, },
        { 0x2e289b1523578dae, 0xf7d9ef7871f6e01f, 0x8d660848663aa29a, 0x0ffad6679dc8931c, },
        { 0x0da98e0693b961f1, 0x12244ea74eab4629, 0xaa09ac0d4b3814de, 0x10fb7f2f160e7061, },
    }, {
        { 0x759da5dfd952abd2, 0x78203568cb26c4fd, 0xfb64f1c68637eff6, 0x5efbe968084e0e5d, },
        { 0xa952791296b7c16f, 0x742ea24463b993f8, 0x7235f33561eae9c1, 0x2e5eac3b7585a4ac, },
        { 0x1b726b3891a749c5, 0xb5586aa3dab40e67, 0x88099153580cc602, 0x3cc1afa4c51bef09, },
    }, {
        { 0x6f976a4aacd6c32e, 0x305329e5ccf43ac4, 0xfe7d4816086f0ca4, 0x007cf519e2e9b1b5, },
        { 0x3ce79c3e68f99ef6, 0x325b0aa5d9e039d6, 0x82cf5b6bfe49c5b4, 0x056e559e1bc2886b, },
        { 0x5a7df532e4fdf1ee, 0x36f0dcbbff4da4c2, 0x61cb3fae1f12517b, 0x5590f3acdeb2ecc2, },
    }, {
        { 0x4172840a806f5fdc, 0x80cf2e21e1c0a075, 0x084594c3491792d5, 0x3563a178f2d73480, },
        { 0x7cbd6b009452d82b, 0x38c36968fd61aff5, 0xf1e876aadcc1b605, 0x2be54637ef44ff57, },
        { 0x2521d70fd94df57a, 0xc116688438867486, 0x35555e4cfe418af3, 0x5e7e4ddbd2b4fe50, },
    }, {
        { 0x863c123742a5e164, 0xf804441040fe626b, 0x9eaca6865ee943c0, 0x440be3ec90eaa131, },
        { 0x5913d7bfcf866787, 0x94dfd2a45537b933, 0x26e14edcfe16e927, 0x15024e35c1e88451, },
        { 0x6e258350c56b7601, 0x01142ceb6aa0350a, 0x3fc64777704ebd76, 0x212d6d6f6c64969a, },
    }, {
        { 0x3f6a655c2270f29e, 0xae2799fc77ab8171, 0xfd162b369f8a1e3a, 0x47ecd96a90a5ce4a, },
        { 0x77a8becc7409dc66, 0x431ddd349d3ed1d2, 0xc62942fa86e266f8, 0x46e5d5ceff103bc2, },
        { 0xef425960a228b0b5, 0xba65476dfcc50da6, 0x41739c393cb98b55, 0x0cad13fa2154e94b, },
    }, {
        { 0x9ee97be636820d13, 0x2c08b8ff70164070, 0x6de8f2e3fa54210f, 0x79107ceab360f28b, },
        { 0x4c23ef35782d32c8, 0x9767dda24a46b57a, 0x72bd6514d4eba29e, 0x299a9aadd2eac7e9, },
        { 0xd01ae2d386f7386d, 0xffd4b450219bb9ea, 0x702d677da8774383, 0x6ff363515329101b, },
    }, {
        { 0x0e0c47457b3587a1, 0xc842158291769720, 0x91003f56fe5627c2, 0x6a2aec01e85ca61c, },
        { 0xfdfd082da1a52002, 0x4ccd6f41eac7ea40, 0x22c95fc6130b9c5a, 0x768ac9e51ffb2792, },
        { 0xedb5a4a20063c3c8, 0xdb4937bdb73fd668, 0xe2f7b92422d13cf7, 0x1778db95503454ec, },
    }, {
        { 0x0f46d1ab7bd82d09, 0x1a9308e39e94d4af, 0x589e4fed09848426, 0x24b29da0750656d6, },
        { 0x9191f66c6d03ffe7, 0x526d9ad8b277c2b7, 0x9198592e0ca63f40, 0x5a85058770681952, },
        { 0x4543aa9d682446b9, 0xa6873bd0d95ee3da, 0x1352667c05e3da81, 0x728a483817331e0e, },
    }, {
        { 0x77ddd129d60a2084, 0x3bc7adee8e3a7568, 0x3d228c358ed6ecd9, 0x54b5c9a949d97f1d, },
        { 0xe8967bfd5108c512, 0x74fb7c212e2a9b05, 0x7c9a0d253a686776, 0x7de8499e0f17bf89, },
        { 0x9fe207c4e7ea37a3, 0x7f881301962c41d7, 0x0db13977582df624, 0x37baa92e3c04302f, },
    }, {
        { 0x6eb6315ad6bb0652, 0x94458890bb775227, 0x98fc0a5b239c7bdc, 0x599baa3ed4ccab9a, },
        { 0x266b206c8adcdb1e, 0x68c72f11fcaa28d2, 0x16eb534e8c0ae8b4, 0x5c3de1c09f0f7437, },
        { 0xdd43e25a33666249, 0xa4f8dd6c8415de01, 0x5503235fabcc2daf, 0x2e43c2f7b80b04f8, },
    }, {
        { 0xa3e3b738140e43e3, 0x51847226c8f5a0db, 0x3ab9407df9ca4102, 0x6db8bb5b83005256, },
        { 0x1d5544ca2b9a38dc, 0x24b58309f0ec30c3, 0xb8882607f41bab99, 0x11dbaf2a48a09817, },
        { 0x1538256f39e7695c, 0xe7e39640b7290d01, 0xaa24a1ead2582546, 0x1b0c351b76a6ab17, },
    }, {
        { 0x7550ec45b78294ad, 0xee25b2b4ffae747f, 0x65613f77c8bd244e, 0x544f0f82ff601899, },
        { 0xc68b1d070c0270a9, 0x6681ae6c986ea7de, 0xe8dcabbd5d0c0161, 0x5aa30e062e487553, },
        { 0x8e41da8d2d9b8e57, 0x67a7574c0728e94d, 0x1b6310d8b69717b4, 0x338cc7be18e8e439, },
    }, {
        { 0x7ef87e29dd00b10a, 0x46d7f7ffa69e8401, 0x10a99a875576be17, 0x2beabc55f16e56b4, },
        { 0x868f3f7f6449ee75, 0xcd76900055b7a8af, 0x4abcca728db11e10, 0x2d193f9f0ac3050c, },
        { 0x7c5ccf4722824c89, 0xaa79b00107baf05f, 0x05e488ee79f0e84c, 0x59cda4e562c70edb, },
    }, {
        { 0x404b413e33da2fde, 0xc3f5a183f66c3990, 0xe33f99a06dfcf1a7, 0x541d736f22cc7500, },
        { 0xc8f2565026bd487a, 0xa80ac2efb1851d00, 0xe30bea78739adca1, 0x13f5094dc6d25f17, },
        { 0x67c7d3374b9ff45f, 0xc49e7a0f5b287bb4, 0xb6f2a7588e115fa2, 0x66e18672da626d0a, },
    }, {
        { 0xbb2e00c9193b877f, 0xece3a890e0dc506b, 0xecf3b7c036de649f, 0x5f46040898de9e1a, },
        { 0x739d8845832fcedb, 0xfa38d6c9ae6bf863, 0x32bc0dcab74ffef7, 0x73937e8814bce45e, },
        { 0xb9037116297bf48d, 0xa9d13b22d4f06834, 0xe19715574696bdc6, 0x2cf8a4e891d5e835, },
    }, {
        { 0x226343e5955ec3d0, 0x7efc48d89d175a1b, 0x98385297a18b1a3c, 0x5492bea1636c4322, },
        { 0x6224329d0d0edba7, 0x467780808e04a966, 0xc62ce016618d99f6, 0x2f810078408f8d31, },
        { 0xe9d5690a9065f588, 0x68c463a035160408, 0x2db1cabb9e8baa87, 0x2e19e362fb39ae79, },
    }, {
        { 0xb131a87740b81edc, 0x2529cce1c50a6e37, 0xf87178c43d86df40, 0x3aa294b972690814, },
        { 0x31746c23083ac1c0, 0x3dac7b511c3db42f, 0xd4f12b16228a9f80, 0x40730ab9b1806c0a, },
        { 0x152402f20552a87d, 0xe69f2c7546352008, 0x09d51a513da75873, 0x4e71cd08428677ab, },
    }, {
        { 0x34ed76d0cd63a5cf, 0x74dfc1951f9660cf, 0x54f23ff7e1127868, 0x09ef2d3f52d2d691, },
        { 0x9473f1cd00004a14, 0x7acdb1fd2d3e101f, 0xbadf67da3f14d3b1, 0x720361d4dd5fa9b5, },
        { 0x75e9d3ef07b12573, 0x07776372d4a699b6, 0x25745c88b4d16ddd, 0x2c42d3a4e7729db5, },
    }, {
        { 0xfb4fb8bf87a8ed1b, 0xe61dff953562c0c1, 0x7e0099fb2eb7c3f6, 0x4bfd51db8338b347, },
        { 0x2cb2fe4fd66daba4, 0x10a1807089aa7f00, 0xb4f37238714aed31, 0x205b2e2c4e9867f4, },
        { 0x3fd78085bdb17451, 0xf4aeaaa7ef06458d, 0xc871854b8c5836f4, 0x0ffd138d1115b50d, },
    }, {
        { 0x9df506730af1f0d8, 0x75535041b11f9b70, 0xe8c3368099ec8b4d, 0x250ad1f281ca5e23, },
        { 0x0580f7ec3e11be9a, 0x0973d59fc90a1d21, 0x9f32c03c33e48155, 0x4aaca15627b19d76, },
        { 0x03d16c211577e699, 0x89a469648e5f2c2c, 0x86374d4cc56666a5, 0x741323b846bb7004, },
    }, {
        { 0xec4f3deda9b25300, 0xcbbaba10d085e701, 0xc24f4daea1518c1f, 0x4dbe8c6ffa6eb46d, },
        { 0x4091c10304a28695, 0xeba17e70b40db506, 0xf5bc279cdd1102cc, 0x1664358723c81ccc, },
        { 0xa24e517143273ad2, 0x60f4788244a356ef, 0x1e1754b3bb30b209, 0x2de4a5e282488c52, },
    }, {
        { 0x2ff8ba76c7f31ce4, 0x71ccb65873599089, 0x5b79d7204b8ff7e3, 0x201e2da6d7e029fc, },
        { 0xbbf5bbf89958d609, 0x30981d35e1a1359d, 0xe326cc62b78d1571, 0x154c89fc4bb2de35, },
        { 0x950868cba485927f, 0x5a65136c99309dbc, 0x2ea9e9f269c7a60e, 0x18f69d0c7431685e, },
    }, {
        { 0x084abc56d78789ac, 0x5659a9504c34a837, 0xfcbc29ccc8e32d7b, 0x490bad9213a35552, },
        { 0x9b3f9942e038d4e9, 0xeda655ec6cd96aa5, 0xb2a94883f0b1d81f, 0x38639ccf3a65a94a, },
        { 0x2af87c80a7817c8b, 0x3290308bcd31c019, 0xc389cfaea7be274b, 0x6a125aeb10bfe889, },
    }, {
        { 0xe20267b63d24e162, 0xca58692468a6d572, 0xd498c9e703adc37a, 0x202d45269bd5eb1b, },
        { 0x67ff338994aaea83, 0xbc661418df086f32, 0x4f6e34a28278692b, 0x3405d3a9171d88ed, },
        { 0x0349f89cf3d71f71, 0x6c59ec91d391896e, 0xc601b73560e3b489, 0x37eefb40d16a68c3, },
    }, {
        { 0xf56e41423ba3c672, 0x573e3974e9368956, 0xd6a2df8ee7a5b7b3, 0x3147f8664a134a97, },
        { 0x961cd0528c5fb918, 0x8120a7a9645ebc80, 0xd448452a0f1e6013, 0x05d53aca06aa1a4c, },
        { 0xc6c34026daed85e1, 0x3cfb0bd6b2d4b1e0, 0x0f12416c370b87cb, 0x0224a2ebaa945a13, },
    }, {
        { 0xe526547cb7294205, 0x0b8e75ed78f2614f, 0xd06a714bf001d89a, 0x6af7c9394efa41c1, },
        { 0x67fabbaefa893cc6, 0x6701db39849aa0e9, 0x01718cf68bc4c43a, 0x0e5f734bdc8b76f5, },
        { 0x0540e82d3cb66a17, 0xd8d957ef8d07d837, 0x9bbd0223d8c809e2, 0x57d04a0ac06c287e, },
    }, {
        { 0x305c39e37e348382, 0xa6961d62151f8ff5, 0x9b3ba96e3e75b0c5, 0x1185b2c041a21da3, },
        { 0xe20ac913080a91b1, 0x1f5176f97f570a48, 0x5ed8ea61add984b8, 0x79613d640fbfcfbb, },
        { 0x7211bde42a195a23, 0xd6b973d161b1b2bd, 0xba6eacdced88d400, 0x04c90947f6cb238d, },
    }, {
        { 0x7f39ca1a829f065e, 0xbe497d51b1aed2ec, 0x7cd059b844cf5f3f, 0x69bccd333df51b03, },
        { 0x210a1fb29bf83ff4, 0x8a6a0f760a04da21, 0xf0c749bf5f656141, 0x00439d12eb8bc7d9, },
        { 0xae9a651b3d9e0c3a, 0x0ab0d68d6d13fd33, 0xfe6951a3035ab3a2, 0x29a6a9cdc39650a1, },
    }, {
        { 0x39fd96e57b6b016f, 0x5abf92ed078b082f, 0xc905f0e6ff012931, 0x762d88a2ee141164, },
        { 0x2a6c73c8c8d474cc, 0xbdd881574240c83c, 0x6482d3d9bfb840ca, 0x42e006c7577fa278, },
        { 0x2bc628bb7b427448, 0xc30a541b2e34ad6d, 0xbf692762d6008dfd, 0x7be4cba5b43988bf, },
    }, {
        { 0x0610b7c9d185cae3, 0x22bcf7f5393e1eea, 0x078f20b650b5a498, 0x528d0c4e75dabca8, },
        { 0x1cbb13cc7b7d5d33, 0x48c2338c897f473a, 0x92d5a13d38dab8c7, 0x73b88c7492780aaa, },
        { 0x0f091f8080e27fb5, 0xd839aa99b49773aa, 0xd83cce6ae18637c9, 0x039465879c92e15b, },
    }, {
        { 0x089f56a4fdcdba71, 0xdd032610923939be, 0x2ffd1803cd3fecde, 0x3508c3552852b233, },
        { 0x976aa0140c6804c2, 0x37c4080304b12d0a, 0x10287e281f0473ff, 0x3d2910dddc487f40, },
        { 0xe0156e4e8c0e2780, 0x96963476147b72ad, 0x35d6e9be15f7d223, 0x7a8ce344732f5bbc, },
    }, {
        { 0x49192b179510c49a, 0x8375825d463d436b, 0x6f7494040bf30fba, 0x32a3650efb3957ba, },
        { 0x875624d2686d3b9f, 0x7d2508ccd010ec63, 0xe63ebb260c70a6e4, 0x5a0317781ae53b8a, },
        { 0xe15d831cd70c0232, 0x21810b5efef230e6, 0x062533b6443cec71, 0x7950df2a3ed57767, },
    }, {
        { 0xe6180436237d9833, 0x525696efeb38f533, 0x09b10fcd5c27bac7, 0x608ea1cf3ba1b69d, },
        { 0xa056ac6fe717a6e0, 0x00a1b3117b5f666d, 0x5e3c449279035a10, 0x72c0e858faf62b07, },
        { 0x07d035bbd207e130, 0xaf7053b550117ede, 0xfeec01c7912c5c81, 0x1a1fac04886e8d51, },
    }, {
        { 0x1912ee7cf01e858e, 0x592c8bcc9acd0fa7, 0x753f805d8450bdcf, 0x71bed11670af3575, },
        { 0x01184df2c2b406d2, 0xf76b3c566f206147, 0x57139091a3822f1e, 0x5cc2b0d9793c5d31, },
        { 0x8f19fbd2d85fd4b2, 0x59f1a77fd4d2ce75, 0x85a6259ec51ef80b, 0x7e12a2de0c738fd2, },
    }, {
        { 0x559ad49cd34f712f, 0x4ea4d1efd8fb35e5, 0x5909ad24f723baac, 0x3183c84ac4704f23, },
        { 0x4ec3d6c0fe1998a0, 0xb9be0fed7fc3bccf, 0x6f7cc86fc92a1f39, 0x0a10c4dffba3b48a, },
        { 0x31fd07b545995411, 0xc9fabf61c069d22d, 0x766b85bbe0252beb, 0x10e588a72188a218, },
    }, {
        { 0x74d498cf82ffef83, 0xece755a2fecd2263, 0xc281538af493ef3e, 0x54d171ebf2303fb5, },
        { 0xa63c9137563dd509, 0xf7efae0715c99b21, 0x64780364a6d3fc6d, 0x1f33151851d6089f, },
        { 0xb871ceb125c78add, 0xc0d7481e3c04c296, 0x7de886113a0f2833, 0x6e1cf0a45b7672f8, },
    }, {
        { 0xe78948aaba0b03cc, 0x46a3c56158b6b3f2, 0x0e9640661be481f6, 0x3664e01c8c25c3ef, },
        { 0x198599bba0a886b6, 0xb6bd8fe7876645f2, 0x2417725c8b7591d0, 0x3a08aae7ee1935d1, },
        { 0xdc1cbf911534af52, 0xf4fc0066f7af099b, 0x97397e9a88137747, 0x349b9669e80e74f4, },
    }, {
        { 0xdd240e0aef471ffb, 0x57702113a7856ef8, 0xb8e69c9fd957ef9c, 0x79add1b2bd3c4500, },
        { 0xd90137dd80847f79, 0x71d9affebaf8bfda, 0x8b749f1e4f245bcd, 0x642f1cf2b2ae8bf6, },
        { 0x9744e5e5706c0e4c, 0xb0e8006d24fe9139, 0xa69baed0c7e8ad1c, 0x4ccb9693d6ffb0d6, },
    }, {
        { 0x7e2177f1245d71bb, 0x4bed3a32e461cd75, 0x48eff412c81e1c67, 0x0759d1eabaf0a0e1, },
        { 0x694f0a87ac5573a1, 0x86f3edac8d490435, 0x8d05aa8c69cc424f, 0x3285f2cdd3c5308b, },
        { 0x615777a6cb0ca64e, 0x4b2bcc26f5a307b6, 0x7e2a4b4e819e835f, 0x60a0d80bd949088f, },
    }, {
        { 0x033e6d4240daab62, 0x07c1b46cd01be421, 0x826ff014ad6293b4, 0x7d7bb82e65a1773a, },
        { 0xbe064448510d12e8, 0x76ae88e19c8cba8a, 0x14357377591be232, 0x28a272d71b5d7b44, },
        { 0xfe3ab0d4869e6f11, 0xa5125ffbfa8b0bfd, 0x6b7075c9453b4759, 0x0f03ff742f32f6a8, },
    }, {
        { 0x97949b173c94139b, 0x24983ad3077870c7, 0xf170d9e9a52b4d16, 0x1427c7b4d359bd78, },
        { 0xc227d0a39de5470b, 0x0dbb91658f7b3712, 0xaf2ce7802c7feff1, 0x60df4ea5ae5534ff, },
        { 0xd153615398ce5195, 0xa48ebef16a72c08f, 0x4f62191039cf6449, 0x0284779d0a9fbf8a, },
    }, {
        { 0x451e38f30497bedf, 0x24cf56e1895e5de6, 0xc2b0a79b92b28e61, 0x47360aee44ed80de, },
        { 0x04defea1d88c8292, 0x839768c3b3cdd5de, 0xe6adbb6552dac506, 0x1b7c8b7e7eac3d52, },
        { 0xb1a2c1439e1212b4, 0x2a2101c119527854, 0x2aac510b604e4c45, 0x61a9cc16fe1736a5, },
    }, {
        { 0x79fd59e3b14f5809, 0x2267267fd56a8126, 0xef28666f2bbb00cb, 0x04e6471770348032, },
        { 0xe8e260c9c1977f4f, 0xa0286c1829041772, 0x38e868bbcb1f6094, 0x5ec36955a701f2d2, },
        { 0x0f66a28baeae57d6, 0x03629de7cef451d2, 0x261b291936d14cda, 0x61079d97dc7d2447, },
    }, {
        { 0xc2b8f3999b8f0086, 0x410ace9f9868469c, 0xcd834c270a12c56e, 0x004d7ba35d423dd2, },
        { 0x2a6177d845813608, 0x09b400f5e38dffef, 0x1754a95f805633d9, 0x514113b957e41e84, },
        { 0xddaf0ea9dd05fa4e, 0xfc381a60b53379a2, 0x2d0f18dfc2d4c1c4, 0x271fa439e1897bf5, },
    }, {
        { 0xe252270fdcebe6ad, 0x1e1d8802049efc21, 0xaefddec98a4eafc2, 0x2f35657f0119a165, },
        { 0x7833fda9b3ffca5d, 0xcd8f493a3c5cd47b, 0x3a68f5f418647670, 0x57dea19ad69044aa, },
        { 0x4d255ced967ed429, 0x50651b4c822a1e3a, 0x9a93f1c0a0f56f58, 0x26da368d2a033c0f, },
    }, {
        { 0x6416f27a4ed2fcbf, 0xe2b4cc903a355dde, 0xd110663f6f7cbc57, 0x46c5572b41497ccb, },
        { 0x17e813d3a7b1d798, 0xf9a26a8f6b76bcc6, 0x124c1bb731e5ca85, 0x71b4d3008bbaebd3, },
        { 0x1fd1c0aa83ede334, 0xfa257c0985cb7b8c, 0x13ccc95b96f02600, 0x5ee86aecc80c67b0, },
    }, {
        { 0x77102749f3fddfe3, 0x093b131779beac52, 0x2c1a5cc1baba54c7, 0x7cc2effd82ed20f9, },
        { 0xfe5e30fe627c7b1c, 0x7d6de30d8d82fc2f, 0x6b7981b11c1bf394, 0x6ef79537955690ef, },
        { 0x317bd4c938ac1d9d, 0x120b22abd43a48c8, 0xd373691a25ba47c7, 0x444f44ebdde79e2d, },
    }, {
        { 0x7ca77213d05f2d3c, 0x8a61a13f0c380077, 0x3b6956cac2c06c3b, 0x21ee01f2aab52726, },
        { 0x7aa57251b1f4422d, 0x210f89842356b185, 0x8db759f1093e1544, 0x0e5204246ccb7827, },
        { 0xe99ccc6c7210d454, 0x9278b2198dbed73a, 0xc404d1c3c984a72b, 0x14afa6b532746a36, },
    }, {
        { 0x6ff950ea2bf4376a, 0x970d30f65b86db4a, 0xfd04f109d2ee3eb1, 0x23654dedb4b8924c, },
        { 0xbb24631daa496a92, 0x9a4dccaa560f5b8e, 0x964c615bfad2d247, 0x0cf07c2d65681c6a, },
        { 0x6b1e6c89d9512bc6, 0xff307335c7992ae0, 0x1e41b43978620616, 0x5f5527f073477c30, },
    }, {
        { 0x0f59397b366f2064, 0xca579e0400ffbdfd, 0x4627c0746dbd3e97, 0x101ae639eb010236, },
        { 0x72743e1adf923a45, 0xa0bd821a1b73554c, 0xa643c9635cb5b0a0, 0x48ea9eaf01173dd7, },
        { 0x42bf06ed1472d9e0, 0x6c5b1064572013ef, 0x71536c8628f3a65c, 0x0aec303795c00785, },
    }, {
        { 0x6cbcddec843606bb, 0x47208c7a1b031e0a, 0x804710d373c8628e, 0x24767bce8c6389bb, },
        { 0x5921911d4deab9f3, 0xb29306f0d2ee2c08, 0xebe05d5236213424, 0x6fc95025487fff46, },
        { 0xeae74b3752e0d95a, 0x88fc646d01319900, 0xa4648c1c6d6482b1, 0x185ef4b1c69f70c4, },
    }, {
        { 0x07bb7b82626330a6, 0x290663cb102512be, 0x3512692369fcc45d, 0x6360f82570f829aa, },
        { 0xda1c29b7673dea49, 0x5c791fa27acdb58d, 0x5e538aac691cd92d, 0x19d69f88986252b8, },
        { 0x82021d93adfc8d6d, 0x9dd0fdb3784c8063, 0x383f70a2e64adcce, 0x5f3913384333e288, },
    }, {
        { 0x5af519089ad6dad0, 0x3f696e0b3f313daa, 0xc73864cf87370a41, 0x729d8f120ef8379f, },
        { 0xfc02938c3c528cb3, 0x8a0c04b4139eac04, 0xd580e5fb6303418a, 0x58a3ef09ac77136c, },
        { 0xb13fa0df25cf504b, 0x787038f05187842b, 0xfc5e7ddff3929b12, 0x4d698742cd706af1, },
    }, {
        { 0xf08f2f01dbe7c0fb, 0xeab65c3e6b294170, 0x761ef854f7d940dc, 0x0ed7ad78ad5649ab, },
        { 0x3c18c4f5a94fabd9, 0x09a0b15c5d4decbd, 0xd93cef9b767627d1, 0x46b61cf57dfa8237, },
        { 0x4324e59d12288fd0, 0x2cb474455329a7b3, 0x3b3ad463d9e94476, 0x47be3633a3243d8a, },
    }, {
        { 0xa38193936bd68f25, 0x921a6473fc692aa7, 0xfbaab570102c552b, 0x22260e700da1e96d, },
        { 0xd2dd61e7dbc441a8, 0xf5a113db3004a0f1, 0x49f05e44f70504ce, 0x2a59cd573b07b9a0, },
        { 0x7d785536e116ab23, 0x2d2ffa28e3fee820, 0xf8d8fff8c62f75ce, 0x7626171b15c34f43, },
    }, {
        { 0xd6887915c7f57505, 0x03f31b6122d43e1a, 0xeba0accd4c4c6736, 0x028d5c21839a3696, },
        { 0x4c2b8a631aa0566b, 0xceeb045f7df7c625, 0x5a5c972e2bda1196, 0x6b8b1fb5631d1ad0, },
        { 0x855229f02d0cd5b7, 0x5965aa97c0948473, 0xd7b3004f7766caff, 0x0d9d4eb5fcbcc3a6, },
    }, {
        { 0x2e9b3e1c2c8cdc4a, 0xf301bddb39292b8f, 0x908063db6f45b448, 0x03d2a9b125c80dff, },
        { 0x3ecadb0a72684fb0, 0x6ddfb9ed18c4f682, 0xa94247e46dbbbaef, 0x2c63df06adb3d562, },
        { 0xfd788e43ba8dd8fc, 0x30fae1c25eb4e217, 0xbd60c01e00568551, 0x066f5248f35c4364, },
    }, {
        { 0x780d2dc1d486eefe, 0xddf38e2964578dcb, 0x1246a3ab7430fbbc, 0x4d86a09ba8cd4997, },
        { 0x9bf5440c101afa1f, 0x2ff63482d49d4402, 0xe31ebee76e5af731, 0x3bbdf993417f08f9, },
        { 0x441ea7b8546b1ac3, 0xf37d7e19483a5d2d, 0x00edf00eba1e5334, 0x2aa79702144b8647, },
    }, {
        { 0x32c35992bd027592, 0x040cf4a5ca53e5d6, 0x89b6026ae0463dd4, 0x52d53fff62f20b12, },
        { 0xac40c7ef8c16314b, 0x9b3d64e9c5723c00, 0xe9d3de862ed1b06d, 0x6f30202862ad53d1, },
        { 0x174eeab8d94095dd, 0x613ff392a5969d4f, 0x045d9a4f8e5a87b0, 0x02b2bd7ea9223a7c, },
    }, {
        { 0x6881a0dfb995530d, 0x6e83a1c84bed5faa, 0xdc3a5d25a2c68798, 0x147ed1907558ce6a, },
        { 0xdb936784762dd969, 0x465ba33cd457ba9c, 0x355491a691eaf814, 0x7b1caaad4cce2b18, },
        { 0x72311a1a04ff82ff, 0xdf23622663cb8d77, 0xb5355304d3f8974a, 0x7411c190f77d8471, },
    }, {
        { 0x6ca9d966064182ac, 0x1f3fb4ae49d1aa8e, 0x52d36be5e1cfa08e, 0x1d5f55868add0c82, },
        { 0xd776b5cabe030055, 0x912f1b63f8c69c37, 0x493ab77cbf914632, 0x748a0934cc3cda8d, },
        { 0xd0a445073e4787e2, 0x59559ab0df698de5, 0x75d355566859c01b, 0x0b3b867f11fff4bb, },
    }, {
        { 0x7ebc968669d85f9f, 0x73f2592293669d5c, 0x0e003301bc8c1531, 0x711add89995c9c5b, },
        { 0xad4a3101aafab7ad, 0xf2957758684bce24, 0xa2d2ca0aa594eef7, 0x1089cbcbefb47df6, },
        { 0x031bf98091dbd7bd, 0x516d988ece5d7df8, 0xab63762d85288be5, 0x7e70e158eeb64b4e, },
    }, {
        { 0x8558038dc22fc777, 0xf1d4eed59d6072de, 0x08290ffbe6005486, 0x4a3a9ac23da24681, },
        { 0xdcc797fd2f1518f2, 0x619d3a39b7963f13, 0xffee1a8cf4a7a3b8, 0x35b6adfeb05fc1bc, },
        { 0x4252edaa94c5a6ab, 0xf43c48baa251c3ed, 0xf62e9831f4ae7037, 0x0f5e3b9db185b22f, },
    }, {
        { 0x04a3b229a77b1206, 0x73fdca1232729905, 0x5c40e6458c540ccb, 0x486fb19627160146, },
        { 0x6e4fa0009f8877d9, 0x86ae024c94dc3bbf, 0xd8f74efb5fbf9c12, 0x7fa9e45bed93e8f0, },
        { 0xaa44b5b79f91b00a, 0x7191aaddc8e0e05d, 0x381527f8209afdee, 0x6a143f29338ccba6, },
    }, {
        { 0x535c23e68efe5f3a, 0xc95059e56118bfad, 0xb4077ebbc89ee4a2, 0x3a5905ce3953b360, },
        { 0x02e8b328b2eb2a77, 0x23b9d06da5a66d2d, 0x0e778bc30f2c344d, 0x5c916ef3a2531b50, },
        { 0xc2edb16896c8547a, 0x37284976f32371b3, 0x87f4fcc0521bd71e, 0x25e95a25518d4d5f, },
    }, {
        { 0xa784cfbdc59cce86, 0xde20df51efe84612, 0xc74dc0c3d6d09035, 0x1ce27d48c62f2772, },
        { 0xdeb818f06900db2d, 0xc3e266cb38dd986a, 0x1e22df7f5b65ff3a, 0x6e3cda36d08609c2, },
        { 0xb67fe90bea68db73, 0x31de909ba3a5b89d, 0x12e3020a3df076f6, 0x6f329a5b417f0ab1, },
    }, {
        { 0xaa18e903eb9e5375, 0xd02c3cc006b0aeba, 0xd3ce459eb5a31ea4, 0x570ede91d1dae60c, },
        { 0xdc239ae1c5c320d1, 0xbc34b40df867c4e4, 0x172d49e5b870c507, 0x033f1c738eb1f1bf, },
        { 0x6d90fbf279d70924, 0x9e2843d34f12e885, 0x84c95865e8f279c9, 0x3acee246b2593b63, },
    }, {
        { 0xa16ccf7c7f00e7eb, 0x44d48a1c5d2f2453, 0xe6a70ebbfe858d77, 0x1c7e8edd00c3f6cc, },
        { 0x5a4124f476f6f17a, 0x5016ea68ca560ccf, 0x74aa0a4fc33b4658, 0x5e5d38adb36d2380, },
        { 0x10e9863f3383a423, 0xee18b967505cfa1d, 0xee2c240290a3dcf0, 0x695425ffa9ffe46e, },
    }, {
        { 0x2ca6c19089dda3dc, 0xcb471b4d41dc9ead, 0x871f7b14dad19a35, 0x68aeed648765ce86, },
        { 0x5d43801f88964a3b, 0x4729a1ac63c01e8e, 0x060b46bc2aba8d3b, 0x1e0df677dd261858, },
        { 0x7bd27536df365e9d, 0xd1fd7dc406b14bf0, 0xee218b17ce4a5e1f, 0x78b98eee8fd166a5, },
    }, {
        { 0xa7f04383e899a9ef, 0xea99f9f269457248, 0xb999d0c0eb3c8884, 0x71d2b0901c6a5651, },
        { 0x323e86345203596a, 0x35149575c2ce6bf9, 0xb34240750a01df70, 0x24ba40cca34c1460, },
        { 0x3a8695c3b864d752, 0xa918a66e8a15f21c, 0xb568e03ad8591bf3, 0x71a0109b349dce00, },
    }, {
        { 0xbd43151079f24806, 0x185162c4631cab26, 0x804d93737f839e2b, 0x6c92e2184aaa443d, },
        { 0x3258fe15cfa352e1, 0x3729a9635da8a720, 0x094e2ee036cf728a, 0x3cce260b92288216, },
        { 0x6bdd19138e60f089, 0xa496bbee78bd3671, 0x51cd50469ad2f4fe, 0x5b4188a3400836cd, },
    }, {
        { 0x2e950a832a0a8ec5, 0xcb8f99321e681a60, 0xdd9d800552279d7b, 0x14a74efb09c83204, },
        { 0xa3733c5562318b07, 0xafd22928e83ea099, 0xeaa02a3bd27d5c8c, 0x33c2a00d235c0e94, },
        { 0xb2946dd44b576dba, 0x21e413dbd8b53f91, 0xca632a5f2be355cc, 0x71ecbe0ecb961edf, },
    }, {
        { 0x407aed1c545707de, 0xd93c6610674278b5, 0x16b3bdde309c6338, 0x41a4ab631125bfbf, },
        { 0xa46dc40b1764bad8, 0x48a207e518d9bdd6, 0xb82ee200ac541818, 0x3a1ff80edf85a2ce, },
        { 0xb21902439c9860ed, 0x1613409ecf7f2882, 0x37220e85c178472d, 0x32cfb2aeac8dbd3a, },
    }, {
        { 0xb25e771029805697, 0xd7cdef18ca8fb69e, 0x04c6dfe4abc42411, 0x778b5db87cc589c1, },
        { 0x33438f9d05afc320, 0x7a5733609fadd4dd, 0xb7d737323f1b9745, 0x6be47b8dc77145bf, },
        { 0x1f178641d61248bc, 0x416d5a3bc6cd5075, 0x879d3b3446ba559e, 0x72ed94850f83e801, },
    }, {
        { 0x1638b3932d59f526, 0x013ea2ecaf1ce3e3, 0x510a1fcc4228d072, 0x388aad93fbb99466, },
        { 0x420ced3ebe203ea7, 0x64542dcee565ebb2, 0x2114d04d0396360b, 0x07c5598f22edac34, },
        { 0xff6534d580586512, 0x4fffa5a4fad70ec0, 0xd4f1bb112b0eb2c0, 0x506ed52f0f3597fb, },
    }, {
        { 0xe534a5a178638c57, 0x30eb2b0a0cd5ca53, 0xed53f11bf9b39956, 0x7c102dc875581672, },
        { 0x003c16fda045f0cf, 0x52ce729fd76cd5ce, 0x3ac6c289458a2035, 0x26584c06845ce344, },
        { 0xe0404454a6fa9a6b, 0xadbd84574792baa1, 0x6620eef32c58377b, 0x37a562bab0896918, },
    }, {
        { 0x557a8b50d44592e0, 0x680e77648142e64f, 0x3d778549c8265d83, 0x028ccac098101691, },
        { 0x0cc1bc852794fca4, 0xb76c830e25a418a0, 0xca12bcf62c7d2445, 0x254dc7532ebfd751, },
        { 0x22fb9b36d6f443bb, 0x0fbd567e83bebc96, 0x07f1da95e71352c2, 0x7af3b15425b2ac21, },
    }, {
        { 0x8d12481a85407e5a, 0xbff0af2e7c7da4f3, 0x4b5294921f78174b, 0x3a341284929caf23, },
        { 0x6dfdf91cb37ea56d, 0x844fe3d5a6322508, 0xcccb3c9ba01c7e92, 0x0a1eeb119ce2ca2d, },
        { 0xd14afb05b8dce9d6, 0x936127702086abcb, 0x28bd9665978a592a, 0x508171f7f04dab40, },
    }, {
        { 0x8ef4f973cd52cc12, 0xbb0df0e302fd31ba, 0x5e14df026ff47a51, 0x09867330be09d912, },
        { 0x471cecd2d63f6ac3, 0x0d189e50e006f7e5, 0xc5d653d77c53490d, 0x11994587720a1d58, },
        { 0x159456a8042830c3, 0x4adde4baec968cee, 0x1eecbe25ea1f266f, 0x28b6fd853cb564c2, },
    }, {
        { 0x83c99a1e0f77f81f, 0x6210474eb468d0e2, 0x47f5c79de6c8c914, 0x1ef81e9f930fbb04, },
        { 0x5468751940d54f15, 0xf25dde44dc0fb02a, 0xef69233af03dd8e3, 0x118f2566a1bb8f0d, },
        { 0x578b7a8d9341b419, 0x6268cbe60957a946, 0x935236cc87aa2f02, 0x2776a859286055d7, },
    }, {
        { 0x48143d7a329b2fdb, 0xbff719e2c223aa1c, 0xb31bb4d29c3a8deb, 0x0f8d065c634c45ad, },
        { 0x4c9f1fc7ccde3be1, 0x2435505590c32eb4, 0xc580cafde2d2f279, 0x6f6c8fafc5e561ec, },
        { 0x8a6ce79ca1722ba6, 0x2e0e52714d71598d, 0xdd784d487cd608ec, 0x6d6a760565dd681f, },
    }, {
        { 0x68d9f28b67d0431b, 0xc76d28149cf3a930, 0x56928f210d71bcb6, 0x171c31ccbeae381e, },
        { 0xf90bcd5c03aafb2f, 0xd8f88081815934c6, 0x1f78f418f0d10cff, 0x02a080ddd6871cac, },
        { 0x8f144aad8e0fa607, 0x6fbfda7b59157879, 0x1a501ef4ad515cec, 0x6b4088ac2280ddb3, },
    }, {
        { 0x969142aa8050a9b3, 0xcefd6d4a57911e94, 0x3d2d1a0775b77969, 0x147f20fbe8d20f48, },
        { 0xca043e96be55ce98, 0xa6faf419b754ae8e, 0x4e7b4f8a184d0c43, 0x1f3e7406627e08b2, },
        { 0x4260f29a6d195f8e, 0x74048da6a5e9ba6c, 0xcdc8f2c125e1b139, 0x038ab2c7629ddbbe, },
    }, {
        { 0x1f79eedc7248e519, 0x6b3156e5021df483, 0x7b41679092559021, 0x31cf66c10c02bb33, },
        { 0x80641606a7df6991, 0x16e4a154983c009b, 0xaef405e033439b23, 0x38dbd5fa3b59d939, },
        { 0xea714c2c4fa75ad0, 0x611ac296ad4bb08f, 0x7fc710eeb94f3a1c, 0x0047a91289a71dac, },
    }, {
        { 0xbc270ad52eeae5b2, 0x960205ab88a30147, 0x83d7459f77fda1a5, 0x2df4e0519a174564, },
        { 0x31235505c6892a5b, 0x7a9bf56bf1dfa44b, 0x741af28fc5166795, 0x683d4e51e0acfef7, },
        { 0x21afa2444d9f5cc8, 0x703c9308fa1e2e70, 0x9f846432c54d90dc, 0x6f0fab9f2671bd0a, },
    }, {
        { 0x9306efdf8b6d507a, 0x8a858e5824a4a4af, 0x8d5239cdc3ae6195, 0x7378ccd4128d97ab, },
        { 0x685b27dcc37e5d11, 0x5a9823afe6053657, 0xa7aee36887c245f4, 0x5f6b1088d4dd09a1, },
        { 0xa06accf6f5bf02f1, 0x0d0771784884ce04, 0xfae1f25b72fb231d, 0x49daf8819d8846e3, },
    }, {
        { 0x474bd6817c4ae902, 0xf002198a43a99bf3, 0x8dfbc6bb64dd2e2b, 0x2254069ec02c6d59, },
        { 0xfdecfcfb6dd1c7df, 0xd91b028b0a18991d, 0x74e4247fc4d0cd0e, 0x2f95c6945b8b5add, },
        { 0xabed41391f7a3922, 0xb6a3683363649dc3, 0xaa21a8b1a91c5e11, 0x2d9cf4386d6027c5, },
    }, {
        { 0x50468b7624388bc3, 0x851ce5fd087e73d2, 0x36ae4762191c2aed, 0x0a60bbe0b011a87a, },
        { 0xad19cfcaf225f04c, 0x845dff33deab77a4, 0x0d9fa0befdcaaf08, 0x5e52efa6841efe92, },
        { 0x4545a8c74c7dd5c8, 0x5bf56fffb14a069e, 0xc76c1af786429545, 0x670c07e93055a690, },
    }, {
        { 0x8f9eb10044b91294, 0x67478749e1542fed, 0xb29b469b535dfb48, 0x264a3066e2e25f76, },
        { 0x0cb27f60f83f9eb8, 0x56adff6c8cd82b0b, 0x037ccda8b03fb08a, 0x42183695d2cf0873, },
        { 0x0b717ac4a872b3a1, 0xa15a62d23d88d1b6, 0x922e432426d780da, 0x732113b26e905fe5, },
    }, {
        { 0x1ff3a45b7ebbfbb9, 0xa8d86bf7d1282547, 0xb0d70c8c1be1d08c, 0x06d2a8cfb21d19ed, },
        { 0x254ab26b340a1440, 0x3b55db81b7faf723, 0xd46c27206f1b0db0, 0x78efb14cd54f357e, },
        { 0x68e9bd0eba181d54, 0x13b4d6bf3d5a5ab0, 0x00b66ea20eae1fae, 0x265d02eedbe5e1bc, },
    }, {
        { 0x2a0df4f8b50c5365, 0x804eb0a8c691a617, 0xd243f6cea48d86d7, 0x16276e2427f6fe1d, },
        { 0x4992e44e7a6c53e4, 0xa53ff2b1ab8b7de2, 0xcd0aaacd165ad0de, 0x30f58c44608625e4, },
        { 0x866d31d0b5136309, 0x7e074b405593d3bd, 0x29c71c281f16aec0, 0x2ebbacae573fb2f2, },
    }, {
        { 0xedc5296f3555106d, 0x4d434b2d8321a7e2, 0xc0070068e894362e, 0x6582a9616d6b4542, },
        { 0x213e29a3fb09c2a0, 0xb609d3797508a69a, 0xc18c2ebb350ae36c, 0x1e187f5b2b4bd901, },
        { 0x6d2241cbb07eb9ca, 0x49644e096f5e89ce, 0xc54d97c8738c0d4c, 0x53fe0d2e4ad173b9, },
    }, {
        { 0x97f3291d71b8eed7, 0x18c9e9914789e136, 0x08550516c8cdaaca, 0x05395ec388719efd, },
        { 0x8d59447d2690ef5f, 0x7bbbe17c8d9006d4, 0xef2330170e7a90bf, 0x2a748b9302f22abe, },
        { 0xbbeaa6fe3402ca3a, 0x5faa5a4297bf31ce, 0x15b0d145029baadc, 0x5b8923891e21f9f4, },
    }, {
        { 0xaf3675229aa4b0d6, 0x20c3be00761ee007, 0x03ef18166338301c, 0x47bab7a2936bd202, },
        { 0x3b1e3a0b7bb8e420, 0xf05b0cbbb63a34c1, 0xaae01ae4f7980bda, 0x64b18fb2409046b7, },
        { 0x0dfbf734432ffe0b, 0x3dbee455ec8ddd34, 0x4a74149f66ef9396, 0x63a250d15ed6e73c, },
    }, {
        { 0x98e51b908128f86b, 0x263f13e5e6dd7e30, 0xd06ace319689b811, 0x573293484f23ee94, },
        { 0x024fcff51a4756d4, 0xb020630a6b882e91, 0xb5d5f7032075f533, 0x324509d078b94de6, },
        { 0x6775cd9962140ea0, 0xd642d377365b3ef0, 0x601901ba3f65cb56, 0x2b939509cf83f05b, },
    }, {
        { 0x7e5e99ff6d6a253c, 0xd97f397a30a1733c, 0x670e3805ab987888, 0x030d83cbf8852227, },
        { 0x8b92fede6e1239e2, 0x71acff21161b082a, 0x237048d284b89c7c, 0x0ae5961865c1d212, },
        { 0xd8daf9cb8fd54819, 0x39e291d5d693e577, 0x182476f6abc8fcc2, 0x4d40af2720b3e13a, },
    }, {
        { 0xa1f94b8bbf688f75, 0x77fb61904405cc0c, 0xaef8d012c85681aa, 0x745c2f628b9e8deb, },
        { 0xdd4f3641264a3ee3, 0x66835f1375fde909, 0xb54a1863d6b2e022, 0x5b427cf0f4ba14fa, },
        { 0x7077ea22a292db19, 0xc78565e544a509e2, 0xba708ad05a325bae, 0x08b85fe992e8d98b, },
    }, {
        { 0x1119ef6d3c2a3fc6, 0xd5ff6c26fd6996a5, 0x77b8cbb20b8a76d4, 0x5effd1c431e72e57, },
        { 0x3c5098ee301d8692, 0x3aaa542c8281b7f5, 0x218c79a2a0c369e2, 0x6419a7ca4c179097, },
        { 0x18ea4f7859c60441, 0x7b83a61699ea645c, 0x8ef0fb5967151f89, 0x5e39c0a8894aace6, },
    }, {
        { 0x0739fb1553bb5de7, 0xa216864abd9b28ae, 0xfb9ad874ba0679e8, 0x2218cba5c38fb07c, },
        { 0xd6ac0ca87cf82652, 0x57869cd56c7fae03, 0xbce832092d89124e, 0x4df95d3fc4c878e1, },
        { 0x59f7733fb5927b72, 0x6be2fea361d984c5, 0x0dfdf83dbae6bfe3, 0x212806be86c56e0f, },
    }, {
        { 0x2e28325d77c9b45b, 0xdf743df4e9846205, 0x85601b41a85f1015, 0x725c2a235e52839d, },
        { 0x087b45fe2640b516, 0x6150594de44c0b2b, 0xba23332a707e28f4, 0x3d262bb0870e6b6b, },
        { 0xa5585ec002400f16, 0xcd8e00eda7fc7c05, 0x5aeaf2b1f56d3635, 0x55defad84d374520, },
    }, {
        { 0x8db456d803fdd76a, 0x5b9e20ffb283770c, 0x26fe8032b98218a6, 0x2f6060ddfb32b99d, },
        { 0xa26cb9e5b5788d49, 0xa71302f4049962c5, 0x2452441f564a2c87, 0x1a7d4aa33fdfc76f, },
        { 0x28072a25009052ed, 0x9a949b5dc3c8eaae, 0x6b8c4c5c1ccfdd1a, 0x7d80ea9fe7d97235, },
    }, {
        { 0x3a8a63cf1ab2966c, 0x213900406b8e081e, 0x75cd30fa0fc10c06, 0x753e0bbebaf97e6e, },
        { 0x8d82476b9b2594be, 0xf40f7b1ea6b5c9c3, 0xe21f8b31a37a28d4, 0x1fb2af1a89c5e37a, },
        { 0x99dbc21eb9d8a7a5, 0xef895a84e18019f6, 0x4e2dfe2c0c75cef7, 0x5be903eb1778ab0f, },
    }, {
        { 0xca7012d21ae40d9f, 0x9871980f30854715, 0xd11ccc26bb143a88, 0x66dff11ff56ec92f, },
        { 0xa10b7043438abf0a, 0xb0d19e865f23e52c, 0xb4199f0e3f32327c, 0x55efb0132c4e9f71, },
        { 0x1b583aa8a0438592, 0x45ca9cdc3ec1bd80, 0x146aa0056283d900, 0x3d4c06cd5fea8e37, },
    }, {
        { 0xa864c021c6427729, 0xd47a8465a8ff6f61, 0x307e2d76bea8cf54, 0x3989530ff1e6aa68, },
        { 0x3a7ca5320c07dcb4, 0x594847eb555ea8bd, 0x49b689089291248d, 0x52542cddf3238c74, },
        { 0xc8a55dbdaa940777, 0x84fc512ef602d143, 0x5837062a7afd5a3b, 0x7c3601ca2c32672c, },
    }, {
        { 0x8eadbc062a98b263, 0x134752d94ecbbea9, 0x3808ef6a674d813a, 0x1439727e8ef4fc28, },
        { 0xf05c508a7a068932, 0xd1e9f2c75073d890, 0x54748a79c75e6dde, 0x307685c74c059128, },
        { 0x4fe76dc1c15dac95, 0xcec2c5348c08202d, 0xf1476e540a41138a, 0x35961d0ee3e1d667, },
    }, {
        { 0x4bbf1d5f9b209918, 0xcc66599f45b9707b, 0x6836ed77e8027a49, 0x14fdf6a658725c1d, },
        { 0x7035b550c5584125, 0x47cadbfda4396714, 0xf5e3d12d842f7962, 0x5bba73f7eae771cb, },
        { 0xbe0fcdf3d88b7241, 0xef43e3ee12c1c616, 0x242ce3379520131e, 0x780900a85142263f, },
    }, {
        { 0xe15b351eeec0b00a, 0x05ebde80f40560e2, 0x54fc15314bf3d79c, 0x25acf601420e3b9c, },
        { 0xf607d3a53b286403, 0xdf4cac2d1bcbb5a3, 0xc3f2619ca728de9e, 0x4e6bc727d4610144, },
        { 0x1c635c8c9fe2d6fb, 0xe7045e088d4c5e93, 0x1fec2422d4dbc6ac, 0x7a66cd4c2cf86cb9, },
    }, {
        { 0x6b03097a3650ee88, 0xd648b90949aa5c27, 0xa71404fb2eacca57, 0x0b3e2f2e1d0b2598, },
        { 0xd5c62b19429b19a0, 0xbe74d35bedb7fbc2, 0xcb5ecd511501b39c, 0x520afa73cd68f30c, },
        { 0x13d10d7d9b94f78a, 0xdd41952cf275fceb, 0xd07381c1faf7428d, 0x5f35320c554e6426, },
    }, {
        { 0x1d5fbb2b42498232, 0xaf964cd764617b01, 0xc73bd158683f21e0, 0x3f5bccc5c6a76a1f, },
        { 0xe869823290e427f3, 0xc6b82be5f566e95f, 0xc5515d073ab95dd5, 0x7fd106bed23b78d7, },
        { 0xbe1af04a3f0206e9, 0x2a63d67881aafadd, 0x28aec8f0c0797d2b, 0x5e60065bc1cdf03a, },
    }, {
        { 0xab40f73de1a07256, 0x530709301dca0bc8, 0xc0cf5ce6cce7012d, 0x382861fa7db88cec, },
        { 0x2628f9cd8f6343a8, 0x9a9d3259ec504b17, 0xeb95f901211e7aa1, 0x2a3c28620ce7a6b4, },
        { 0x1fc22389b50d49a2, 0x016c7980598ce6ee, 0x1d36f19a766bda55, 0x0ea73a031a3397cd, },
    }, {
        { 0x3d34e768921f8aab, 0x3746bea36eee4795, 0xee02d1aa59fa3cf7, 0x3273e234daef126b, },
        { 0x5740ca72d35a83cd, 0x0cde016001abea1b, 0xfb7a7f6aecc21b16, 0x01bf2d13b048ac23, },
        { 0x3c4a8bba71ce9ab2, 0xf73c2eb3b61c7208, 0x73859ef5b8f9e211, 0x6a321b746d310cf1, },
    }, {
        { 0x604e1e0467109e02, 0x8085e5ed08f2d4b9, 0xe0b5d551c9753421, 0x5f89495a3fa55db6, },
        { 0x75d6c8b487e6e9ca, 0x843f057b036b42c6, 0xb93110bc0cbbb0cf, 0x142a78f8ec836a5d, },
        { 0x5092782d9efe3033, 0xd22832a0989329d1, 0xa47c589130f98d17, 0x0c8a16518b27a411, },
    }, {
        { 0x8c7bae652d02e7bc, 0xea3fa9ee0cf6b274, 0x2380ce11cb7cf44d, 0x0d4b1e8e2cae4701, },
        { 0x6bac21a99e42a45a, 0x30e796db951eafa0, 0xfb9ba2d32d65766a, 0x61c09de270457fa6, },
        { 0xc8a7d67a0c9487b9, 0x485058e7ffde5271, 0x8d4c7c77900e3d73, 0x65fbb4fad77ce7fd, },
    }, {
        { 0x2133f962a7a629e8, 0x5c9cb54979a0d1b5, 0x73904bb0065f8c08, 0x4c4ed976fad0e729, },
        { 0x51a01fa8aa47117e, 0xb2ebbdb0024b909c, 0xea98d10f5a972df4, 0x67c48354e0a11b75, },
        { 0x3a993910dbf660ef, 0x45e54784c8da689f, 0x7bb3de28bd91ef05, 0x402ff4ecccf3b029, },
    }, {
        { 0x1d35d625d6230ebb, 0xd633770352db2b63, 0xeb98655849d27e43, 0x6dfe341f222e4495, },
        { 0xc69f489b08b7bdbc, 0x7781295483679268, 0x8e7db0bac79b057c, 0x516f4adf38daafe8, },
        { 0x1252ddf7063d2b1e, 0x3bbac44decc83029, 0x35b88786d75b2461, 0x67a510b7bc3fcd97, },
    }, {
        { 0x72b274bc76e5a91f, 0x740896fde4e2fb58, 0x3db49cbc4398d057, 0x27877b9767bf7c11, },
        { 0x2aa0fceb249ad55e, 0xa6f928ed1eddd923, 0xc8e4240243a7f8a7, 0x615d6eeb5403a69a, },
        { 0x1ad54899e331cef9, 0xcbc382f26b042736, 0xc1fa864cd766bbaa, 0x1bd72747d8807d0b, },
    }, {
        { 0x5a0ee465d41947d6, 0xb9280af4d666a96d, 0xd511a413481e9302, 0x1f446d472085229c, },
        { 0xf3f6bdf787182e8c, 0x1425737b4ef9320d, 0x6c45a7be82810bb9, 0x3abe9d28a27d07a6, },
        { 0xa797b1f73dc9c9d1, 0x8779591b6557305d, 0x5d50c61286ca651f, 0x25a59857570c11c1, },
    }, {
        { 0x30ce8ac0efa598fb, 0xd4bc57be9fb323f4, 0xcd8ab81c25e14200, 0x51311426d81fe2ac, },
        { 0xf43ff7b28eb09162, 0xfa5b38c34b73fd5e, 0x6a23d0b89017c42b, 0x455bfc850cfd6207, },
        { 0x63b74d5867d08476, 0x3132dc5759b9f88b, 0x01aeac8cd1c53732, 0x4f799cb2fb5116b1, },
    }, {
        { 0xc22d01280a37b21a, 0x49ed37dd70d97402, 0x3d37f157d66b2847, 0x6a2b356a0aec3f5b, },
        { 0x39680d7771af1bbb, 0xba9eb095466dba37, 0x381319b3d8e2a8ad, 0x64e6a9ab21e7d92d, },
        { 0xe63dc0eaa2cc3d15, 0xb16beda73cbf1fa6, 0x1c68e13cf673e230, 0x7966a776416aba4c, },
    }, {
        { 0xb558aba43d82da00, 0x47d3988b1ab1358b, 0xa19837407929eba2, 0x2fbadd1fd0ba5596, },
        { 0x787d73033c7abd11, 0x69196d20cf018842, 0x6f0507167fbafa24, 0x066548ff83651bae, },
        { 0x07bc0bbfce0f5195, 0x3b4a0144a136c3b4, 0x13fdef15877423bc, 0x15b8d97ba9a62435, },
    }, {
        { 0x10842629439b69a5, 0x23eaf65014f24042, 0x9587f46c6e0b55df, 0x7cc8e2711e1543f4, },
        { 0xb6b9da7674252ca8, 0x5299ec18a9c3ed11, 0xf6a7059d0fd15b94, 0x0b39d3428e8de417, },
        { 0xcf729d4e7356040b, 0xa187a711b3a71f12, 0x017fbccb876626fd, 0x4c62d119314c6e1f, },
    }, {
        { 0xf6e6c5e0929993e0, 0x848d3e2ed0acf6a2, 0x9cadd0487be75b7e, 0x05d76fdb9872704b, },
        { 0xa067f291415f92e0, 0x8213233b866a63d5, 0xa310bd28dbb61c6b, 0x58c7f6584d19886e, },
        { 0x0d9c95c818432a70, 0x7b43333f6fbded5b, 0x1a86acb8032f2ba8, 0x1420e549b3c5e5c4, },
    }, {
        { 0x6462d6a753dce318, 0x24d923d6f7b16df5, 0xad2e2e7602a99281, 0x2e6dba46277fbe41, },
        { 0xf9b31ceff3cf928f, 0xe330d91afd306987, 0x2c212c31bc042d05, 0x27414b23414d6a2e, },
        { 0x80c26d33413b240e, 0xa53193aea0b89246, 0x4e1a8487b2eaf649, 0x12082216e3b4c292, },
    }, {
        { 0xff3f7f45f2e69ef1, 0x5504f08fe2ad789c, 0xb0fb7e50aa245200, 0x166f05bdc6fe94e0, },
        { 0x1f29fac801805a82, 0xa4fea30766a06cfd, 0x05f50ee61e96f1b1, 0x43a107e2e4685e28, },
        { 0x0597cd661d28d94c, 0xdc3b2d1e9ec745f2, 0xf8257aa6e546e123, 0x02d21113117c9a30, },
    }, {
        { 0x848082598f37293a, 0x0b0a87fadb95552b, 0x73963f882e59a403, 0x36d75cdd8ebd7260, },
        { 0x4d242941c258f5d4, 0x8dd203e74015ed2c, 0x008c9429949c2c89, 0x0052b3e338ed2a3c, },
        { 0x7aaaf491abfe286e, 0x97d15582a4ec546a, 0x5d4eb2951e9810ef, 0x46f21171e8784f3c, },
    }, {
        { 0xd6e2272ef130cd06, 0xbe6b24d6a5fca011, 0xb243eb1f0cc83c25, 0x3b6640990bc8fc42, },
        { 0xc009f143900fb0c1, 0xb3bdeec454ecfbaf, 0x22ee7d00db64d788, 0x026b096233ce3df1, },
        { 0x12f86676202dccf2, 0x732a8168b996e643, 0xf54ff9dc4c422c52, 0x26b11500b75c7a10, },
    }, {
        { 0x26788efaf191a8dc, 0xe8dfccfce7263590, 0x157c9362aa2026d1, 0x4e32caad4a5d144f, },
        { 0xec9b891d563adb50, 0x8d78b1db1f671dce, 0x69617115282d197c, 0x1ae2643a5573a978, },
        { 0x09636fb74fcf2434, 0xdd881e20af995a40, 0xce211a895f91c2b7, 0x6a50999cfd4900aa, },
    }, {
        { 0xa901a64e916d6651, 0xb968b5a2c8cf8f39, 0x5d33a4a9b8aa4e43, 0x2f5ea7466e1e792d, },
        { 0x89a32e63e8d4a1b7, 0x81f8e9db4d6cf0c6, 0x216cc7b74039b74a, 0x3b19c5629abf60c3, },
        { 0xe81e14498737934e, 0x6b3499b665c399f4, 0xd7978199b51e6884, 0x105bc430b4895825, },
    }, {
        { 0xe96b4aa0f621f64f, 0xe10d4effe5323071, 0xcfbcc47d14199139, 0x05f714e8e03e31c0, },
        { 0xe990577a12a2ca2b, 0xc6bfa1a781c27939, 0x9bff68dcad3f21b2, 0x2ed1957140615194, },
        { 0xb0f981e77feb4a3d, 0xdd112bd4137c3028, 0xf8f3ebc588da7346, 0x380c6fbb676e052a, },
    }, {
        { 0xe23740a63f0d7841, 0x4f44e3a0d710a7fc, 0xb0e2cdb55c773e57, 0x2767944e242057a1, },
        { 0x1788cbb3dccf6fa1, 0x60664abe5ca89daa, 0xf2615418b2d8d3f5, 0x2cb8cc249a87a9aa, },
        { 0x76f79f1b49a9a2e4, 0x7761490d0d3e6934, 0x94eca394e338299e, 0x4c66d193423fa351, },
    }, {
        { 0x298607417b167fb3, 0xd8ee7462a96f4d50, 0x915860d005fa8a5c, 0x7a2dbab46abfe43e, },
        { 0x055e414815756360, 0xb92d55def13335e8, 0xdf4f2f240e7a7e61, 0x6ec868ca40d9cb0d, },
        { 0xafe1d51d99cb6197, 0x6cd0be1308d20839, 0xe193e85d57381600, 0x4a122e9be7bf8571, },
    }, {
        { 0xf5acb918e945b433, 0x1ae37ffa79fc5fe1, 0x911bc333ed7ba805, 0x6b8f213afae89e01, },
        { 0xa5ab485be3f61f56, 0x6bd0e91ac5bcc9af, 0xb149a0bce1f420cb, 0x29486d5d047fc073, },
        { 0xb2d4da5ad4d2f3c5, 0xd0e3f287e50034ab, 0x1aec0bb185ef0211, 0x7a17675cfce4de1e, },
    }, {
        { 0xb96b2855706f8116, 0x1382f1d67ccc4854, 0x532583ad92a438eb, 0x7166252117757420, },
        { 0xf0d8fb4fb79c5ece, 0x85ed103e7a504844, 0xf0d85b2fe20450c9, 0x44b642cfd72f1427, },
        { 0xa55dae4de06c0ae0, 0xd7d3ef74931cdd9e, 0x8bb2e9e87eb1eb19, 0x1f23796302376355, },
    }, {
        { 0x4f29137a2e03c97e, 0x6f7520cc3bc118f1, 0x3c3678a9007645d9, 0x3f939e7f759754ef, },
        { 0x0b4a4b89f3408157, 0xd91b244d98bb6e01, 0x3c9cf5cd92f96920, 0x055c9eb8ff26e2b2, },
        { 0x555539a7321fc3bf, 0xa30e2c0585be27f6, 0xe7cfae16b912f285, 0x5cc126cc88c26918, },
    }, {
        { 0x0e8307eb90d994db, 0xe7bb920a6e890359, 0x16af416a35e90551, 0x3a28c2ae13531c29, },
        { 0x677bcef0e1054bc8, 0xe84d0129c1a3cc47, 0x7b1c8dd2736efa12, 0x4c9f0b71553fa1b2, },
        { 0xf928da27302a6af3, 0x45821985765eed6e, 0x1da41c9a9f16381f, 0x7f7571495a8b665e, },
    }, {
        { 0x338c379b6ecea774, 0xd70dadd800e29b08, 0x63f15548424479ce, 0x79116d861250212a, },
        { 0x4bd56a4f01d4a737, 0xd5f275e740aa9569, 0x5164aefe03e1a57a, 0x6a7834c3243bd46b, },
        { 0xb11b32f25bba7838, 0x8e04c3bed630ddc6, 0xa8fff1716fcc10b6, 0x17002bd08b186b7d, },
    }, {
        { 0xa3ab22984ec4be53, 0x6576e77e2ca3ce50, 0x741be0f6e33f81ff, 0x222cac2ca2b99ec0, },
        { 0x1bf02d8cd9d254d9, 0x060693fc0fd963df, 0x6e475f0fdfcf5d5f, 0x3a4ccec9a566b430, },
        { 0x5972d89428233251, 0xeac181198be4a49e, 0x22a092970d68bf4c, 0x0d6013f7feaa3b93, },
    }, {
        { 0x2cd6a70ceb928119, 0x2692ede445bf3f92, 0x8c5caca790bf49da, 0x7e4d74b6e24118d8, },
        { 0xa18731c1809dcc58, 0x11434d054030ef96, 0xe68da9829fcbd010, 0x636ee3facc0c179a, },
        { 0xcd9ce414c3e9bb85, 0xc10e619eb637120a, 0x75be15219a75ee40, 0x30508bfa2735afef, },
    }, {
        { 0x4d53005316bb0d37, 0xbac24d5203cc5483, 0x1f28e20c5dc0cb4e, 0x1630777e8d0672e0, },
        { 0x95f5141366e30df1, 0x5f411ec09a0fdc58, 0x20f50861563ea93a, 0x43003be0af225922, },
        { 0xf2273b68bb5a5dd1, 0xdaf11c932811bbe7, 0x953a0a86f51b2b74, 0x09ff24c1c7580cc5, },
    }, {
        { 0x82dc02b14f3bb0ba, 0xd9941c8fc31e249e, 0x2c7a4850632d7957, 0x35228a4f7599b474, },
        { 0x3e20bab847ec2a0e, 0x40022ac38d305391, 0x2a093a82c65117fa, 0x4978381c3fafa224, },
        { 0x59ac5693f6b009c7, 0x38ebc33eb989b083, 0xbfaf9a7cbfa62030, 0x2e233e3a3bf74ab3, },
    }, {
        { 0x667d59d91f661a2f, 0x90c091b0e0b13add, 0x4509a47fb9cfeecd, 0x390e4209ebdc5481, },
        { 0xc672031460dbdd3d, 0xfd6bc073c22c8799, 0x46fd56cb94f8c734, 0x79ab3492ee0f1343, },
        { 0x0b2b4597b8c2a508, 0x71f8fbe0ba5a1c4f, 0x403b3743349ee41b, 0x05041b7a6d7db9df, },
    }, {
        { 0x9e2b2efb94bb6b90, 0xf5d4c66e4a6a3fe1, 0x82af5828d578b4ed, 0x6793d1afa6068797, },
        { 0x9a55c9c68af07244, 0x867f3f9841dd90a4, 0xd4894f77f66845fc, 0x2d01687d414914d7, },
        { 0x61f36185d61358bc, 0xe51bd50fdd3dd058, 0x13d477006c1a8b59, 0x16c7949a98b5a8bc, },
    },
};

static const u64_t const_inv_d[4] = {
    0x25e0f276cdc9f843, 0x0b5dd6984279542e,
    0x2b162114cdb9cf66, 0x40907ed214d5ce43
};

static const u64_t const_random_z[4] = {
    0x459e8ba873ff0c54, 0xdbf5f46a84f9550b,
    0x65f44ffd6db0e537, 0x2cf92b822aed4e78
};

void curve_times_base_point(point_affine_t* q, u64_t s[4]) {
    const point_pc_affine_t* p;
    point_proj_t Q;
    u8_t cut[32];
    size_t i;

    curve_fold8(cut, s);

    p = const_base_fold + cut[0];
    prime_sub(Q.x, p->ypx, p->ymx);
    prime_add(Q.y, p->ypx, p->ymx);
    prime_mul(Q.t, p->t2d, const_inv_d);

    prime_add(Q.z, const_random_z, const_random_z);
    prime_mul(Q.x, Q.x, const_random_z);
    prime_mul(Q.y, Q.y, const_random_z);
    prime_mul(Q.t, Q.t, const_random_z);

    i = 1; do {
        point_dbl(&Q);
        point_add(&Q, const_base_fold + cut[i]);
    } while (i++ < 31);

    prime_inv(Q.z, Q.z);
    prime_mul(q->x, q->x, Q.z); prime_mod(q->x);
    prime_mul(q->y, q->y, Q.z); prime_mod(q->y);
}

void eddsa_ed25519_key_generation(
    /* outputs */
    void* private_key, /* 32 bytes */
    void* public_key,  /* 32 bytes */
    /*  inputs */
    const void* random /* 32 bytes */
) {
    point_affine_t q;
    shahf_ctx_t hfctx;
    u8_t hash[64];
    u64_t t[4];

    shahf_begin(&hfctx);
    shahf_update(&hfctx, random, 32);
    shahf_end(&hfctx, hash);

    hash[ 0] &= 0xf8;
    hash[31] &= 0x7f;
    hash[31] |= 0x40;
    memcpy(t, hash, 32);

    memcpy(private_key, random, 32);
    curve_times_base_point(&q, t);
    point_pack(public_key, &q);
}
