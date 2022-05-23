/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq.h>
#include <lcrq_pvt.h>

void test_parms(size_t F)
{
	rq_t *rq = rq_init(F, 1024);

	/* 4.4.1.2.  Source Block and Sub-Block Partitioning */
	part_t src_blocks = rq_partition(rq->Kt, rq->Z);
	part_t sub_blocks = rq_partition(rq->T/rq->Al, rq->N);

	size_t KL = src_blocks.IL;
	size_t KS = src_blocks.IS;
	size_t ZL = src_blocks.JL;
	size_t ZS = src_blocks.JS;

	size_t TL = sub_blocks.IL;
	size_t TS = sub_blocks.IS;
	size_t NL = sub_blocks.JL;
	size_t NS = sub_blocks.JS;

	test_assert(rq->N == NL + NS, "N (%zu) == NL (%zu) + NS (%zu)", rq->N, NL, NS);

	test_assert(rq->Z, "Z(%u) > 0", rq->Z);

	test_assert(CEIL(CEIL(rq->F,rq->T),rq->Z) <= 56403, "ceil(ceil(F/T)/Z) <= K'_max");

	/* verify 1 <= K <= 56403 (K'_max) */
	test_assert(1 <= rq->K && rq->K <= 56403, "1 <= K (%lu) <= 56403", rq->K);

	test_log("F   = %zu\n", rq->F);
	test_log("WS  = %zu\n", rq->WS);
	test_log("Al  = %u\n", rq->Al);
	test_log("T   = P' = %u\n", rq->T);
	test_log("SSS = %u\n", rq->SSS);
	test_log("SS  = %u\n", rq->SS);
	test_log("Nmax = %u\n", rq->Nmax);

	/* calculated */
	test_log("Kt  = %zu\n", rq->Kt);
	test_log("kl  = %zu\n", rq->kl);
	test_log("Z   = %u\n", rq->Z);
	test_log("N   = %u\n", rq->N);
	test_log("KL  = %u\n", KL);
	test_log("KS  = %u\n", KS);
	test_log("ZL  = %u\n", ZL);
	test_log("ZS  = %u\n", ZS);
	test_log("TL  = %u\n", TL);
	test_log("TS  = %u\n", TS);
	test_log("NL  = %u\n", NL);
	test_log("NS  = %u\n", NS);
	test_log("K   = %u\n", rq->K);
	test_log("K'  = %u\n", rq->KP);
	test_log("L   = %u\n", rq->L); /* always 107 ? */
	test_log("P   = %u\n", rq->P);
	test_log("P1  = %u\n", rq->P1);
	test_log("U   = %u\n", rq->U);
	test_log("J   = %u\n", rq->J);
	test_log("S   = %u\n", rq->S);
	test_log("H   = %u\n", rq->H);
	test_log("W   = %u\n", rq->W);
	test_log("B   = %u\n", rq->B);

	test_assert(rq->L == rq->KP + rq->S + rq->H, "L = K' + S + H");
	test_assert(rq->P == rq->L - rq->W, "P (%u) = L(%u) - W(%u)", rq->P, rq->L, rq->W);
	test_assert(rq->B <= rq->L, "B (%u) <= L (%u)", rq->B, rq->L);
	test_assert(rq->L > rq->KP, "L > K'");

	assert(rq->N == NL + NS);
	assert(rq->F <= 946270874880);
	assert(rq->T % rq->Al == 0); /* T must be a multiple of Al */
	assert(rq->kl);
	assert(1 <= rq->K);
	assert(rq->K <= 56403);

	/* the object MUST be partitioned into Z = ZL + ZS contiguous source blocks */
	assert(rq->Z == ZL + ZS);

	rq_free(rq);
}

int main(void)
{
	loginit();

	test_name("Calculate and verify encoding parameters");

	/* RFC restricts F to 40 bits = 946270874880 bytes */
	for (size_t F = 1; F <= 946270874880; F *= 2 + 1) {
		test_log("--------- TESTING PARMS F = %lu\n", F);
		test_parms(F);
	}

	test_rusage();

	return fails;
}
