/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022-2024 Brett Sheffield <bacs@librecast.net> */

#include <matrix.h>
#include <gf256.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <unistd.h>

#define VSZ 256

/*
 * Method adapted from the technique described in:
 * J. S. Plank and K. M. Greenan and E. L. Miller (2013)
 * "Screaming Fast Galois Field Arithmetic Using Intel SIMD Instructions"
 * http://web.eecs.utk.edu/~jplank/plank/papers/FAST-2013-GF.html
 */
static __m128i mul_128(__m128i A, uint8_t y)
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

void matrix_row_mul(matrix_t *m, const int row, const int off, const uint8_t y)
{
	uint8_t *d = matrix_ptr_row(m, row) + off;
	const int max = m->cols - off;
	const int mod = max % 16;
	const int maxv = max - mod;
	int j;
	for (j = 0; j < maxv; j += 16) {
		__m128i D = _mm_loadu_si128((const __m128i_u *)&d[j]);
		D = mul_128(D, y);
		_mm_storeu_si128((__m128i*)&d[j], D);
	}
	for (; j < max; j++) {
		d[j] = GF256MUL(d[j], y);
	}
}

void matrix_row_mul_byrow(matrix_t *m, const int rdst, const int off, const int rsrc, const uint8_t y)
{
	assert(y);
	uint8_t *d = matrix_ptr_row(m, rdst) + off;
	uint8_t *s = matrix_ptr_row(m, rsrc) + off;
	const int max = m->cols - off;
	const int mod = max % 16;
	const int maxv = max - mod;
	int i;
	for (i = 0; i < maxv; i += 16) {
		__m128i S = _mm_loadu_si128((const __m128i_u *)&s[i]);
		__m128i D = _mm_loadu_si128((const __m128i_u *)&d[i]);
		S = mul_128(S, y);
		D = _mm_xor_si128(D, S);
		_mm_storeu_si128((__m128i*)&d[i], D);
	}
	for (; i < max; i++) {
		d[i] ^= GF256MUL(s[i], y);
	}
}
