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

	for (int x = 0; x < (1 << 31); x *= 2) {
		tup = rq_tuple(rq, x);
		test_assert(tup.a >= 1 && tup.a < rq->W,
			"a is a positive integer between 1 and W-1 inclusive");
		test_assert(tup.b < rq->W,
			"b is a non-negative integer between 0 and W-1 inclusive");
		test_assert(tup.d1 == 2 || tup.d1 == 3,
			"d1 is a positive integer that has value either 2 or 3");
		test_assert(tup.a1 >= 1 && tup.a1 < rq->P1,
			"a1 is a positive integer between 1 and P1-1 inclusive");
		test_assert(tup.b1 < rq->P1,
			"b1 is a non-negative integer between 0 and P1-1 inclusive");
		x++;
	}

	rq_free(rq);

	return fails;
}
