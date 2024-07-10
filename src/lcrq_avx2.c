/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022-2024 Brett Sheffield <bacs@librecast.net> */

#include <lcrq_pvt.h>
#include <gf256.h>
#include <string.h>


/* just sum the elements to get r, they are 1 or 0 except for HDPC */
int count_r_avx2(uint8_t *p, int len)
{
	int c = 0;
	/* FIXME - AVX2 should be faster. It isn't. */
#if 0
	for (int vlen = len / 32; vlen; vlen--, p += 32) {
		__m256i v = _mm256_loadu_si256((const __m256i_u *)p);
		__m256i cmp = _mm256_cmpeq_epi8(v, _mm256_setzero_si256());
		uint16_t bitmask = ~_mm256_movemask_epi8(cmp);
		c +=  __builtin_popcount(bitmask);
	}
	len &= 0xff;
#endif
	for (int vlen = len / 16; vlen; vlen--, p += 16) {
		__m128i v = _mm_loadu_si128((const __m128i_u *)p);
		__m128i cmp = _mm_cmpeq_epi8(v, _mm_setzero_si128());
		uint16_t bitmask = ~_mm_movemask_epi8(cmp);
		c +=  __builtin_popcount(bitmask);
	}
	len &= 0x0f;
	for (; len; len--, p++) c += *p;
	return c;
}
