/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2023 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq.h>
#include <errno.h>
#include <inttypes.h>

int main(void)
{
	rq_t *rq;
	rq_oti_t oti;
	rq_scheme_t scheme;
	uint64_t F;
	uint16_t T;
	int rc;

	loginit();
	test_name("3.3 FEC Object Transmission Information (OTI)");

	F = (uint64_t)test_randomnumber(UINT32_MAX);
	T = test_randomnumber(UINT16_MAX / RQ_AL) * RQ_AL;
	rq = rq_init(F, T);

	rc = rq_oti(rq, &oti, &scheme);
	test_assert(rc == 0, "rq_oti returned %i", rc);

	test_assert(F == rq_oti_F(oti), "rq_oti_F() returns F");
	test_assert(T == rq_oti_T(oti), "rq_oti_T() returns T");
	test_assert(rq_Z(rq) == rq_oti_Z(scheme), "rq_oti_Z() returns Z");
	test_assert(rq_N(rq) == rq_oti_N(scheme), "rq_oti_N() returns N");
	test_assert(rq_Al(rq) == rq_oti_Al(scheme), "rq_oti_Al() returns Al");

	rq_free(rq);

	return fails;
}
