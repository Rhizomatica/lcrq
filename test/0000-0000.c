/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <errno.h>

static void test_parms(size_t F)
{
	rq_t *rq = rq_init(F, 1024);

	/* 4.4.1.2.  Source Block and Sub-Block Partitioning */
	part_t src_blocks = rq_partition(rq->Kt, rq->Z);
	part_t sub_blocks = rq_partition(rq->T/rq->Al, rq->N);

	size_t ZL = src_blocks.JL;
	size_t ZS = src_blocks.JS;
	size_t NL = sub_blocks.JL;
	size_t NS = sub_blocks.JS;

	test_assert(rq->N == NL + NS, "N (%zu) == NL (%zu) + NS (%zu)", rq->N, NL, NS);

	test_assert(rq->Z, "Z(%u) > 0", rq->Z);

	test_assert(CEIL(CEIL(rq->F,rq->T),rq->Z) <= 56403, "ceil(ceil(F/T)/Z) <= K'_max");

	/* verify 1 <= K <= 56403 (K'_max) */
	test_assert(1 <= rq->K && rq->K <= 56403, "1 <= K (%lu) <= 56403", rq->K);
#ifndef NDEBUG
	rq_dump(rq, stderr);
#endif

	test_assert(rq->L == rq->KP + rq->S + rq->H, "L = K' + S + H");
	test_assert(rq->P == rq->L - rq->W, "P (%u) = L(%u) - W(%u)", rq->P, rq->L, rq->W);
	test_assert(rq->B <= rq->L, "B (%u) <= L (%u)", rq->B, rq->L);
	test_assert(rq->L > rq->KP, "L > K'");

	test_assert(rq->N == NL + NS, "N == NL + NS");
	test_assert(rq->F <= 946270874880, "F <= 946270874880");
	test_assert(rq->T % rq->Al == 0, "T must be a multiple of Al");
	test_assert(rq->kl, "kl = %u", rq->kl);
	test_assert(1 <= rq->K, "1 <= K");
	test_assert(rq->K <= 56403, "K <= 56403 (K' max)");

	/* the object MUST be partitioned into Z = ZL + ZS contiguous source blocks */
	test_assert(rq->Z == ZL + ZS, "Z == ZL + ZS");

	rq_free(rq);
}

int main(void)
{
	test_name("4.3 Calculate & verify encoding parameters");

	/* RFC restricts F to 40 bits = 946270874880 bytes */
	for (uint64_t F = 1; F <= 946270874880; F *= 2 + 1) {
		test_log("--------- TESTING PARMS F = %lu\n", F);
		test_parms(F);
	}

	return test_status;
}
