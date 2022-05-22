/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq_pvt.h>

int main(void)
{
	rq_t *rq = NULL;
	int deg, v;

	loginit();
	test_name("5.3.5.2 Degree Generator");

	/* hand-calculate some values to verify */

	rq = rq_init(345893, 1024);

	v = 0;
	deg = rq_deg(rq, v);
	test_assert(deg == 1, "rq_deg(%i) == %i", v, deg);

	v = 1;
	deg = rq_deg(rq, v);
	test_assert(deg == 1, "rq_deg(%i) == %i", v, deg);

	v = 6400;
	deg = rq_deg(rq, v);
	test_assert(deg == 2, "rq_deg(%i) == %i", v, deg);

	v = 983913;
	deg = rq_deg(rq, v);
	test_assert(deg == 15, "rq_deg(%i) == %i", v, deg);

	v = 1048576;
	deg = rq_deg(rq, v);
	test_assert(deg == 31, "rq_deg(%i) == %i", v, deg);

	rq_free(rq);

	return fails;
}
