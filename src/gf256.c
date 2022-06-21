/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <gf256.h>

uint8_t GF256LR[256][2][16];

uint8_t gf256_add(const uint8_t a, const uint8_t b)
{
	return a ^ b;
}

uint8_t gf256_exp(const int e)
{
	assert(e < (int)sizeof(OCT_EXP));
	return OCT_EXP[e];
}

uint8_t gf256_log(const uint8_t v)
{
	assert(v != 0);
	return OCT_LOG[v];
}

uint8_t gf256_inv(const uint8_t v)
{
	assert(v != 0);
	return OCT_EXP[255 - OCT_LOG[v]];
}

uint8_t gf256_mul(const uint8_t a, const uint8_t b)
{
	if (a == 0 || b == 0) return 0;
	return gf256_exp((int)gf256_log(a) + (int)gf256_log(b));
}

uint8_t gf256_div(const uint8_t u, const uint8_t v)
{
	if (u == 0) return 0;
	return OCT_EXP[OCT_LOG[u] - OCT_LOG[v] + 255];
}

__m128i gf256_mul_128(__m128i A, uint8_t y)
{
	__m128i table1 = _mm_loadu_si128((const __m128i_u *)GF256LR[y][0]);
	__m128i table2 = _mm_loadu_si128((const __m128i_u *)GF256LR[y][1]);
	__m128i mask1 = _mm_set1_epi8((uint8_t)0x0f);
	__m128i mask2 = _mm_set1_epi8((uint8_t)0xf0);
	__m128i l, h;

	l = _mm_and_si128(A, mask1);
	l = _mm_shuffle_epi8(table1, l);
	h = _mm_and_si128(A, mask2);
	h = _mm_srli_epi64(h, 4);
	h = _mm_shuffle_epi8(table2, h);
	return _mm_xor_si128(h, l);
}

/* allocate and precompute tables */
void gf256_init(void)
{
	/* Table 0: product of y with all four-bit words */
	/* Table 1: product of y with 8-bit words with last four bits zero (i << 4) */
	for (int y = 0; y < 256; y++) {
		for (int i = 0; i < 16; i++) {
			GF256LR[y][0][i] = GF256MUL(y, i);
			GF256LR[y][1][i] = GF256MUL(y, (i << 4));
		}
	}
}
