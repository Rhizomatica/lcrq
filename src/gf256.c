/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022-2024 Brett Sheffield <bacs@librecast.net> */

#include <gf256.h>

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

#if USE_NATIVE
# include "gf256_ssse3.c"
#endif
