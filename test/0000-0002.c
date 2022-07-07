/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <lcrq_pvt.h>

/* RFC6330 has a particular algorithm for generating pseudo-random num-nums:

	x0 = (y + i) mod 2^^8
	x1 = (floor(y / 2^^8) + i) mod 2^^8
	x2 = (floor(y / 2^^16) + i) mod 2^^8
	x3 = (floor(y / 2^^24) + i) mod 2^^8

Then
	Rand[y, i, m] = (V0[x0] ^ V1[x1] ^ V2[x2] ^ V3[x3]) % m
*/

void test_rand(size_t y, uint8_t i, size_t m, size_t rman)
{
	size_t rgen;
	rgen = rq_rand(y,i,m);
	test_assert(rman == rgen, "%zu == %zu", rman, rgen);
}

int main(void)
{
	size_t y, m;
	uint8_t i;
	uint8_t x0, x1, x2, x3;

	loginit();
	test_name("5.3.5.1 Random Number Generator");

	/* hand-calculate some values to verify */

	test_rand(1,1,1, (2375747048 % 1));

	y = 67866723; i = 43; m = 799999;
	x0 = (y + i) % (1 << 8);
	x1 = ((y >> 8) + i) % (1 << 8);
	x2 = ((y >> 16) + i) % (1 << 8);
	x3 = ((y >> 24) + i) % (1 << 8);
	test_assert_s(x0 == 142);
	test_assert_s(x1 == 187);
	test_assert_s(x2 == 54);
	test_assert_s(x3 == 47);
	test_assert_s(V0[142] == 3569308693);
	test_assert_s(V1[187] == 1853869307);
	test_assert_s(V2[54]  == 2253762642);
	test_assert_s(V3[47]  == 2900516158);
	test_assert_s(2423747970 == (V0[142] ^ V1[187] ^ V2[54] ^ V3[47]));
	test_assert_s(2423747970 % 799999 == 550999);
	test_rand(y, i, m, 550999);

	y = 723874628; i = 199; m = 111111;
	x0 = (y + i) % (1 << 8);
	x1 = ((y >> 8) + i) % (1 << 8);
	x2 = ((y >> 16) + i) % (1 << 8);
	x3 = ((y >> 24) + i) % (1 << 8);
	test_assert_s(x0 == 11);
	test_assert_s(x1 == 58);
	test_assert_s(x2 == 236);
	test_assert_s(x3 == 242);
	test_assert_s(V0[x0] == 1843948209);
	test_assert_s(V1[x1] == 2402488407);
	test_assert_s(V2[x2] == 3308488877);
	test_assert_s(V3[x3] == 1522273219);
	test_assert_s(2102720904 == (V0[x0] ^ V1[x1] ^ V2[x2] ^ V3[x3]));
	test_assert_s(2102720904 % 111111 == 56340);
	test_rand(y, i, m, 56340);

	return fails;
}
