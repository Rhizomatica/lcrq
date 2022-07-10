/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <gf256.h>
#include <matrix.h>

int main(void)
{
	loginit();
#ifdef INTEL_SSE3
	test_name("Galois Field 256 SIMD (SSE3) multiplication");

	gf256_init();

	/* verify table 1 for y == 7*/
	test_assert(0x2d == GF256LR[7][0][0x0f], "GF(256) 7 x 0xf");
	test_assert(0x2a == GF256LR[7][0][0x0e], "GF(256) 7 x 0xe");
	test_assert(0x23 == GF256LR[7][0][0x0d], "GF(256) 7 x 0xd");
	test_assert(0x24 == GF256LR[7][0][0x0c], "GF(256) 7 x 0xc");
	test_assert(0x31 == GF256LR[7][0][0x0b], "GF(256) 7 x 0xb");
	test_assert(0x36 == GF256LR[7][0][0x0a], "GF(256) 7 x 0xa");
	test_assert(0x3f == GF256LR[7][0][0x09], "GF(256) 7 x 0x9");
	test_assert(0x38 == GF256LR[7][0][0x08], "GF(256) 7 x 0x8");
	test_assert(0x15 == GF256LR[7][0][0x07], "GF(256) 7 x 0x7");
	test_assert(0x12 == GF256LR[7][0][0x06], "GF(256) 7 x 0x6");
	test_assert(0x1b == GF256LR[7][0][0x05], "GF(256) 7 x 0x5");
	test_assert(0x1c == GF256LR[7][0][0x04], "GF(256) 7 x 0x4");
	test_assert(0x09 == GF256LR[7][0][0x03], "GF(256) 7 x 0x3");
	test_assert(0x0e == GF256LR[7][0][0x02], "GF(256) 7 x 0x2");
	test_assert(0x07 == GF256LR[7][0][0x01], "GF(256) 7 x 0x1");
	test_assert(0x00 == GF256LR[7][0][0x00], "GF(256) 7 x 0x0");

	/* verify table 2 for y == 7*/
	test_assert(0xea == GF256LR[7][1][0x0f], "GF(256) 7 x 0xf0");
	test_assert(0x9a == GF256LR[7][1][0x0e], "GF(256) 7 x 0xe0");
	test_assert(0x0a == GF256LR[7][1][0x0d], "GF(256) 7 x 0xd0");
	test_assert(0x7a == GF256LR[7][1][0x0c], "GF(256) 7 x 0xc0");
	test_assert(0x37 == GF256LR[7][1][0x0b], "GF(256) 7 x 0xb0");
	test_assert(0x47 == GF256LR[7][1][0x0a], "GF(256) 7 x 0xa0");
	test_assert(0xd7 == GF256LR[7][1][0x09], "GF(256) 7 x 0x90");
	test_assert(0xa7 == GF256LR[7][1][0x08], "GF(256) 7 x 0x80");
	test_assert(0x4d == GF256LR[7][1][0x07], "GF(256) 7 x 0x70");
	test_assert(0x3d == GF256LR[7][1][0x06], "GF(256) 7 x 0x60");
	test_assert(0xad == GF256LR[7][1][0x05], "GF(256) 7 x 0x50");
	test_assert(0xdd == GF256LR[7][1][0x04], "GF(256) 7 x 0x40");
	test_assert(0x90 == GF256LR[7][1][0x03], "GF(256) 7 x 0x30");
	test_assert(0xe0 == GF256LR[7][1][0x02], "GF(256) 7 x 0x20");
	test_assert(0x70 == GF256LR[7][1][0x01], "GF(256) 7 x 0x10");
	test_assert(0x00 == GF256LR[7][1][0x00], "GF(256) 7 x 0x00");

	__m128i A, yA;
	uint8_t out[16] = {0};
	uint8_t in[16] = {0};
	const uint8_t y = test_randomnumber(0xff);

	/* first, lets check our tables are loaded correctly */
	A = _mm_loadu_si128((const __m128i_u *)GF256LR[y][0]);
	_mm_storeu_si128((__m128i*)out, A);
	for (int i = 0; i < 16; i++) {
		test_assert(out[i] == GF256LR[y][0][i], "verify table 1");
	}
	A = _mm_loadu_si128((const __m128i_u *)GF256LR[y][1]);
	_mm_storeu_si128((__m128i*)out, A);
	for (int i = 0; i < 16; i++) {
		test_assert(out[i] == GF256LR[y][1][i], "verify table 2");
	}

	/* fill A with random bytes */
	for (int i = 0; i < 16; i++) in[i] = test_randomnumber(0xff);

	/* load A into SIMD registers */
	A = _mm_loadu_si128((const __m128i_u *)in);

	/* multiply yA */
	yA = gf256_mul_128(A, y);

	/* extract values from yA and verify */
	_mm_storeu_si128((__m128i*)out, yA);
	for (int i = 0; i < 16; i++) {
		const uint8_t v = GF256MUL(in[i], y);
		test_assert(out[i] == v, "%02x * %02x == %02x (expected %02x)", in[i], y, out[i], v);
	}
#else
	return test_skip("Galois Field 256 SIMD (SSE) multiplication (requires SSE3)");
#endif

	return fails;
}
