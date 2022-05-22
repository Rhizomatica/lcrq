/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq_pvt.h>

int main(void)
{
	rq_t *rq = NULL;
	rq_tuple_t tup = {0};

	loginit();
	test_name("5.3.5.4 Tuple Generator");

	rq = rq_init(345893, 1024);

	/* just run through the range and make sure nothing breaks */
	for (size_t x = 0; x < UINT32_MAX; x *= 2) {
		tup = rq_tuple(rq, x);
		assert(tup.a >= 0 && tup.a <= UINT32_MAX);
		assert(tup.b >= 0 && tup.b <= UINT32_MAX);
		assert(tup.d1 >= 0 && tup.d1 <= UINT32_MAX);
		assert(tup.a1 >= 0 && tup.a1 <= UINT32_MAX);
		assert(tup.b1 >= 0 && tup.b1 <= UINT32_MAX);
		x++;
	}

	rq_free(rq);

	return fails;
}
