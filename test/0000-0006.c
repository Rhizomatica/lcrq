/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>
#include <sys/param.h>

#define MAX_SRCOBJ 1024 * 1024 * 1024

/* generate source object of data of random size and content up to max bytes */
uint32_t generate_source_object(uint32_t max, unsigned char **block)
{
	uint32_t sz = randombytes_uniform(max);

	*block = malloc(sz);
	assert(*block);
	randombytes_buf(*block, sz);

	return sz;
}

int main(void)
{
	rq_t *rq;
	unsigned char *srcobj = NULL;
	unsigned char *srcblk;
	uint8_t *intsym;
	size_t F, len;

	loginit();
	test_name("5.3.3.4 Intermediate Symbol Generation");

	/* generate random source block for test */
	F = generate_source_object(MAX_SRCOBJ, &srcobj);
	test_log("block of %u bytes generated\n", F);

	rq = rq_init(F, 1024);

	rq_dump(rq, stderr);


	srcblk = srcobj;
	len = F;

	size_t SBN;
	size_t blklen;

	test_log("ZL = %zu\n", rq->src_part.JL);
	test_log("ZS = %zu\n", rq->src_part.JS);

	/* encode long blocks */
	for (SBN = 0; SBN < rq->src_part.JL; SBN++) {
		blklen = rq->src_part.IL * rq->T;

		/* allocate memory for intermediate symbols */
		intsym = rq_intermediate_symbols_alloc(rq); assert(intsym);

		rq_intermediate_symbols(rq, srcblk, blklen, intsym);
		test_log("encoding %zu bytes\n", blklen);
		srcblk += blklen;
		len -= blklen;
		free(intsym);
	}
	/* encode short blocks */
	for (; SBN < rq->src_part.JS; SBN++) {
		blklen = rq->src_part.IS * rq->T;

		/* allocate memory for intermediate symbols */
		intsym = rq_intermediate_symbols_alloc(rq); assert(intsym);

		rq_intermediate_symbols(rq, srcblk, blklen, intsym);
		test_log("%zu: encoding %zu bytes\n", blklen);
		srcblk += blklen;
		len -= blklen;
		free(intsym);
	}

	rq_free(rq);
	free(srcobj);

	//test_rusage();

	return fails;
}
