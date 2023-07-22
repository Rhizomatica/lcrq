/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

/* params - example program to create RaptorQ context and print derived parameters */

#include <lcrq.h>
#include <stdio.h>
#include <stdlib.h>

int usage(const char *progname, int rc)
{
	fprintf(stderr, "usage: `%s F T`  (F = object size, T = symbol size)\n", progname);
	return rc;
}

int main(int argc, char *argv[])
{
	rq_t *rq;
	uint64_t F;
	uint16_t T, K, KP, N, Z;
	uint8_t Al;

	if (argc != 3) return usage(argv[0], EXIT_FAILURE);
	F = atoll(argv[1]);
	T = atoll(argv[2]);

	/* initialize RaptorQ context */
	rq = rq_init(F, T);

	/* Query parameters */
	K = rq_K(rq);
	KP = rq_KP(rq);
	N = rq_N(rq);
	Z = rq_Z(rq);

	/* free RaptorQ context */
	rq_free(rq);

	printf("F  %12lu\tsize of object to encode (40 bits)\n", F);
	printf("T  %12u\tsymbol (payload) size\n", T);
	printf("K  %12u\tnumber of original symbols\n", K);
	printf("K' %12u\tnumber of symbols, including padding symbols\n", KP);
	printf("N  %12u\tnumber of sub-blocks\n", N);
	printf("Z  %12u\tnumber of blocks\n", Z);

	return 0;
}
