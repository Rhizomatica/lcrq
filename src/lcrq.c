/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <lcrq_pvt.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>

int isprime(int n)
{
	if (n <= 1) return 0;
	if (n % 2 == 0) return 0;
	for (int i = 3; i < n; i++) {
		if (n % i == 0) return 0;
	}
	return 1;
}


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

/* The degree generator Deg[v] is defined as follows, where v is a non-
   negative integer that is less than 2^^20 = 1048576.  Given v, find
   index d in Table 1 such that f[d-1] <= v < f[d], and set Deg[v] =
   min(d, W-2).  Recall that W is derived from K' as described in
   Section 5.3.3.3. */
int rq_deg(rq_t *rq, int v)
{
	assert(v >= 0);
	assert(v <= (1 << 20));

	int d;

	for (d = 1; d <= 30; d++) {
		if (DEG[d-1] <= v && v < DEG[d]) break;
	}
	return MIN(d, rq->W-2);
}

size_t rq_rand(const size_t y, const uint8_t i, const size_t m)
{
	const uint8_t x0 = (y + i) % (1 << 8);
	const uint8_t x1 = ((y >> 8) + i) % (1 << 8);
	const uint8_t x2 = ((y >> 16) + i) % (1 << 8);
	const uint8_t x3 = ((y >> 24) + i) % (1 << 8);
	assert(m); /* must be positive */
	return ((V0[x0] ^ V1[x1] ^ V2[x2] ^ V3[x3]) % m);
}

rq_tuple_t rq_tuple(rq_t *rq, size_t X)
{
	rq_tuple_t tup = {0};
	size_t A = 53591 + rq->J * 997;
	if (A % 2 == 0) A++;
	size_t B = 10267 * (rq->J + 1);
	size_t y = (B + X * A) & 0xffff;
	size_t v = rq_rand(y, 0, 1048576);
	size_t d = rq_deg(rq, v);
	tup.a = 1 + rq_rand(y, 1, rq->W-1);
	tup.b = rq_rand(y, 2, rq->W);
	tup.d1 = (d < 4) ? 2 + rq_rand(X, 3, 2) : 2;
	tup.a1 = 1 + rq_rand(X, 4, rq->P1 - 1);
	tup.b1 = rq_rand(X, 5, rq->P1);
	return tup;
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
	rq->J = T2[rq->KP].j;
	rq->H = T2[rq->KP].h;
	rq->S = T2[rq->KP].s;
	rq->W = T2[rq->KP].w;
	rq->L = rq->KP + rq->S + rq->H;
	rq->P = rq->L - rq->W;

	/* P1 denotes the smallest prime number greater than or equal to P */
	rq->P1 = rq->P;
	while (!isprime(rq->P1)) rq->P1++;

	/* N is the minimum n=1, ..., Nmax such that ceil(Kt/Z) <= KL(n) */
	for (rq->N = 1; rq->N <= rq->Nmax; rq->N++) {
		if (CEIL(rq->Kt,rq->Z) <= KL(rq->WS, rq->Al, rq->T, rq->N)) break;
	}

	return rq;
}
