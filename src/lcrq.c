/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <lcrq_pvt.h>
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

static uint16_t K_padded(uint16_t K)
{
	assert(K <= KPAD_MAX);
	for (int i = 0; i < T2LEN; i++) {
		if (T2[i].k >= K) return T2[i].k;
	}
	return KPAD_MAX;
}

/* KL(n) is the maximum K' value in Table 2 in Section 5.6 such
	that K' <= WS/(Al*(ceil(T/(Al*n))))
	NB: this will always return K'_max = 56403 unless we have
	a small working memory (WS) */
uint64_t KL(uint64_t WS, uint16_t Al, uint16_t T, uint16_t n)
{
	uint64_t v;
	for (int i = T2LEN - 1; i >= 0 ; i--) {
		v = WS/(Al*(CEIL(T,(Al*n))));
		assert(v);
		if (T2[i].k <= v) return T2[i].k;
	}
	return 0;
}

/* The function Partition[I,J] derives parameters for partitioning a
   block of size I into J approximately equal-sized blocks.  More
   specifically, it partitions I into JL blocks of length IL and JS
   blocks of length IS.  The output of Partition[I, J] is the sequence
   (IL, IS, JL, JS), where IL = ceil(I/J), IS = floor(I/J), JL = I - IS
   * J, and JS = J - JL.  */
part_t rq_partition(size_t I, uint16_t J)
{
	part_t p = {0};
	assert(J > 0);
	p.IL = CEIL(I,J);
	p.IS = FLOOR(I,J);
	p.JL = I - p.IS * J;
	p.JS = J - p.JL;
	return p;
}

void rq_free(rq_t *rq)
{
	free(rq);
}

rq_t *rq_init(size_t F, uint16_t T)
{
	rq_t *rq = malloc(sizeof(rq_t));
	memset(rq, 0, sizeof(rq_t));

	rq->F = F;
	/* TODO what is an appropriate size for WS ? */
	rq->WS = 1073741824; /* 1GiB */
	rq->Al = 4;
	rq->T = T;
	rq->SSS = T;
	rq->SS = rq->SSS / rq->Al;
	rq->Nmax = rq->T/(rq->SS*rq->Al);
	rq->Kt = CEIL(rq->F,rq->T);
	rq->kl = KL(rq->WS, rq->Al, rq->T, rq->Nmax);
	rq->Z = CEIL(rq->Kt,rq->kl);
	rq->K = CEIL(rq->kl, rq->T);
	rq->KP = K_padded(rq->K);

	/* N is the minimum n=1, ..., Nmax such that ceil(Kt/Z) <= KL(n) */
	for (rq->N = 1; rq->N <= rq->Nmax; rq->N++) {
		if (CEIL(rq->Kt,rq->Z) <= KL(rq->WS, rq->Al, rq->T, rq->N)) break;
	}

	return rq;
}
