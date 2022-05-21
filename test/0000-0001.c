/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>

#define MAX_SRCOBJ 1073741824 /* 1GiB */

/* generate source object of data of random size and content up to max bytes */
uint32_t generate_source_object(uint32_t max, unsigned char **block)
{
	uint32_t sz = randombytes_uniform(max);

	*block = malloc(sz);
	assert(*block);
	randombytes_buf(*block, sz);

	return sz;
}

void test_parms(size_t F)
{
	rq_t *rq = rq_init(F, 1024);
	rq_free(rq);
}

int main(void)
{
	unsigned char *srcobj = NULL;
	size_t F;

	loginit();
	test_name("Encoding Step One: Padding");

	/* generate random source block for test */
	F = generate_source_object(MAX_SRCOBJ, &srcobj);
	test_log("block of %u bytes generated\n", F);

	//test_rusage();

	free(srcobj);

	return fails;
}
