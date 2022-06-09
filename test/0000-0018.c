/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>
#include <sys/param.h>

static const size_t MAX_PAYLOAD = 1024; /* MAX_PAYLOAD must be at least Al=4 bytes */
static const size_t MAX_SRCOBJ = 328183; /* OK */
//static const size_t MAX_SRCOBJ = 322731; /* FAIL */
//static const size_t MAX_SRCOBJ = 328182; /* OK */
//static const size_t MAX_SRCOBJ = 328181; /* OK */
//static const size_t MAX_SRCOBJ = 328180; /* OK */
//static const size_t MAX_SRCOBJ = 326000; /* FAIL */
//static const size_t MAX_SRCOBJ = 328000; /* OK */
//static const size_t MAX_SRCOBJ = 327000; /* FAIL */
//static const size_t MAX_SRCOBJ = 327500; /* FAIL */
//static const size_t MAX_SRCOBJ = 327750; /* OK */
//static const size_t MAX_SRCOBJ = 327625; /* FAIL */
//static const size_t MAX_SRCOBJ = 327700; /* OK */
//static const size_t MAX_SRCOBJ = 327650; /* FAIL */
//static const size_t MAX_SRCOBJ = 327675; /* FAIL */
//static const size_t MAX_SRCOBJ = 327685; /* OK */
//static const size_t MAX_SRCOBJ = 327680; /* FAIL */
//static const size_t MAX_SRCOBJ = 327681; /* OK */
//static const size_t MAX_SRCOBJ = 32769; /* OK */
//static const size_t MAX_SRCOBJ = 32768; /* OK */
//static const size_t MAX_SRCOBJ = 32767; /* OK */
//static const size_t MAX_SRCOBJ = 32766; /* OK */
//static const size_t MAX_SRCOBJ = 32765; /* OK */
//static const size_t MAX_SRCOBJ = 32764; /* OK */
//static const size_t MAX_SRCOBJ = 16000; /* FAIL */
//static const size_t MAX_SRCOBJ = 15999; /* FAIL */
//static const size_t MAX_SRCOBJ = 16382; /* FAIL */
//static const size_t MAX_SRCOBJ = 16384; /* FAIL */
//static const size_t MAX_SRCOBJ = 16385; /* OK */
//static const size_t MAX_SRCOBJ = 16386; /* OK */
//static const size_t MAX_SRCOBJ = 4; /* FAIL */
//static const size_t MAX_SRCOBJ = 8; /* FAIL */
//static const size_t MAX_SRCOBJ = 1024; /* FAIL */
//static const size_t MAX_SRCOBJ = 1024 * 1024 * 2 + 42;

/* generate source object of data of random size and content up to max bytes */
unsigned char *generate_source_object(size_t max, size_t *F)
{
	unsigned char *block;
	size_t sz = 0;
	//size_t sz = randombytes_uniform(max);
	if (!sz) sz = max;
	block = malloc(sz);
	assert(block);
	memset(block, 0, sz);
	randombytes_buf(block, sz);
	*F = sz;
	return block;
}

int main(void)
{
	rq_t *rq = {0};
	//uint8_t *encobj = NULL;
	unsigned char *data = NULL;
	size_t F = 0; /* F requires 40 bits */
	int rc;

	loginit();
	test_name("Object Encoder");

	/* generate random source block for test */
	data = generate_source_object(MAX_SRCOBJ, &F);
	test_log("block of %u bytes generated\n", F);

	/* API starts here - - - - - - - - - - - - - - - - - - - - */
	rq_sym_t sym = {0};
	rq_state_t state = {0};
	rq = rq_init(F, MAX_PAYLOAD);
	rq_dump(rq, stderr);

	rc = rq_encode_data(rq, data, F);
	test_assert(rc == 0, "rq_encode_data OK");
#if 0
	size_t olen;
	//sym.sym = malloc(rq->T);
	uint8_t *src = data;
	uint32_t symcount = 0;
	uint32_t symbols = rq->src_part.IL * rq->src_part.JL + \
			   rq->src_part.IS * rq->src_part.JS;
	olen = F;

	rq_dump(rq, stderr);

	rq_state_init(rq, &state, RQ_SOURCE);

	/* verify the source symbols match */
	while (rq_symbol_next(&state, &sym)) {
		size_t off = MIN(rq->T, olen);
		test_assert(!memcmp(sym.sym, src, off),
				"compare SBN=%u, ESI=%u", sym.SBN, sym.ESI);
		olen -= off;
		src += off;
		symcount++;
	}
	test_assert(symcount == symbols, "symbols = %u/%u", symcount, symbols);

	test_assert(rq_symbol_next(&state, &sym) == NULL,
		"no more source symbols");

	rq_state_free(&state);
#endif

	/* generate some source and repair symbols and test decoding */

	// FIXME - rename rq_decode_block() to rq_decode_block_inplace()

	const uint8_t overhead = 10;
	const uint32_t nesi = rq->K + overhead;
	uint32_t ESI[nesi];
	uint8_t *enc = calloc(1, nesi * rq->T);
	//uint8_t *dec = malloc(F);
	uint8_t *dec = calloc(1, rq->K * rq->T);
	uint8_t *ptr = enc;

	memset(ESI, 0, nesi * sizeof(uint32_t));
	assert(F);
	test_assert(memcmp(data, dec, F), "verify source and dec dont match before decoding");

	memset(&sym, 0, sizeof sym);
	rq_state_init(rq, &state, RQ_REPAIR);
	for (uint32_t i = 0; i < rq->K && rq_symbol_next(&state, &sym); i++) {
		memcpy(ptr, sym.sym, rq->T);
		ESI[i] = sym.ESI; /* record "received" ESIs */
		fprintf(stderr, "recording ESI %u\n", ESI[i]);
		ptr += rq->T;
	}
	rq_state_free(&state);

	/* initialize fresh rq context to ensure no stale data */
	rq_t *rq2 = rq_init(F, MAX_PAYLOAD);
	rc = rq_decode_block_f(rq2, dec, enc, ESI, nesi);

	fprintf(stderr, "dec[0]:\n");
	rq_dump_symbol(rq2, dec, stderr);
	fprintf(stderr, "data[0]:\n");
	rq_dump_symbol(rq2, data, stderr);

	test_assert(rc == 0, "decoding ok");
	test_assert(!memcmp(data, dec, F), "decoded data matches source");

	free(enc);
	free(dec);
	rq_free(rq);
	rq_free(rq2);

	/* API ends here - - - - - - - - - - - - - - - - - - - - */

	free(data);

	test_log("test done\n");

	return fails;
}
