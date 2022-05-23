/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>

#define MAX_SRCOBJ 4099

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
	size_t F;

	loginit();
	test_name("5.3.3.4 Intermediate Symbol Generation");

	/* generate random source block for test */
	F = generate_source_object(MAX_SRCOBJ, &srcobj);
	test_log("block of %u bytes generated\n", F);

	rq = rq_init(F, 1024);

	/* allocate memory for intermediate symbols */
	intsym = rq_intermediate_symbols_alloc(rq); assert(intsym);

	srcblk = srcobj;
	for (int SBN = 0; SBN < rq->Z; SBN++) {
		srcblk += rq->T;
		rq_intermediate_symbols(rq, srcblk, intsym);
	}

	rq_free(rq);
	free(intsym);
	free(srcobj);

	//test_rusage();

	return fails;
}
