/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>
#include <sys/param.h>

static const size_t MAX_PAYLOAD = 4; /* MAX_PAYLOAD must be at least Al=4 bytes */
static const size_t MAX_SRCOBJ = MAX_PAYLOAD * 10 + 1;
static const uint32_t nrep = 21;

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
	//for (unsigned char i = 0; i < MAX_SRCOBJ; i++) block[i] = i + 42;
	*F = sz;
	return block;
}

void decoder_tests(rq_t *rq, uint8_t *src, uint8_t *encsym, uint8_t *repsym)
{
	const size_t maplen = howmany(rq->K, CHAR_BIT);
	unsigned char symmap[maplen];
	unsigned char repmap[maplen];
	uint32_t rESI[rq->K + RQ_OVERHEAD];
	rq_blkmap_t sym = { .map = symmap, .len = maplen, .p = encsym };
	rq_blkmap_t rep = { .map = repmap, .len = maplen, .p = repsym, .ESI = rESI };
	const size_t blocksz = rq->T * rq->K;
	uint8_t *dec = calloc(rq->K, rq->T); assert(dec);
	int rc;

	/* set repair ESIs */
	for (uint32_t esi = 0; esi < nrep; esi++) rESI[esi] = esi + rq->K;

	/* First, tell the decoder we have no symbols => FAIL */
	/* clear bitmap => we have no symbols */
	memset(symmap, 0, maplen);
	rc = rq_decode_block(rq, &sym, NULL);
	test_assert(rc == -1, "not enough symbols");

	/* tell encoder we have K-1 symbols => FAIL */
	memset(symmap, 0, maplen);
	for (int i = 0; i < rq->K - 1; i++) setbit(symmap, i);
	rc = rq_decode_block(rq, &sym, NULL);
	test_assert(rc == -1, "K-1 => still not enough symbols");

	/* Now lets tell the decoder we have all source symbols => SUCCESS */
	memset(symmap, ~0, maplen);
	rc = rq_decode_block(rq, &sym, NULL);
	test_assert(rc == 0, "all source symbols received, no decoding required");
	test_assert(!memcmp(src, encsym, blocksz), "decoded data matches source");

	/* Now let's make the decoder actually decode something */
	memset(symmap, 0, maplen);  /* no original source symbols */
	memset(repmap, ~0, maplen); /* full block of repair symbols */
	memset(dec, 0, blocksz);
#if 0
	/* copy in one symbol at a time, skipping first symbol */
	for (uint32_t i = 2; i < rq->K; i++) {
		const size_t off = rq->T * i;
		setbit(symmap, i);
		memcpy(dec + off, src + off, rq->T);
	}
#endif
	sym.p = dec;
	rc = rq_decode_block(rq, &sym, &rep);
	test_assert(rc == 0, "enough source symbols received");
	test_assert(!memcmp(src, dec, blocksz), "decoded data matches source");

	free(dec);
}

int main(void)
{
	rq_t *rq = {0};
	matrix_t A = {0};
	matrix_t D = {0};
	matrix_t A_dup = {0};
	matrix_t C = {0};
	unsigned char *srcobj = NULL;
	unsigned char *padblk = NULL;
	unsigned char *srcblk = NULL;
	size_t F, len; /* F requires 40 bits */
	uint8_t SBN;
	size_t blklen = 0;

	loginit();
	test_name("5.4.2 decoding");

	/* generate random source block for test */
	srcobj = generate_source_object(MAX_SRCOBJ, &F);
	test_log("block of %u bytes generated\n", F);

	rq = rq_init(F, MAX_PAYLOAD);
	rq_dump(rq, stderr);

	srcblk = srcobj;
	len = F;

	test_log("ZL = %zu\n", rq->src_part.JL);
	test_log("ZS = %zu\n", rq->src_part.JS);

	/* encode blocks */
	for (SBN = 0; SBN < rq->Z; SBN++) {
		if (SBN < rq->src_part.JL) {
			/* long block */
			rq->K = rq->src_part.IL;
			blklen = rq->src_part.IL * rq->T;
			test_log("SBN %u: long block %zu bytes\n", SBN, blklen);
		}
		else {
			/* short block */
			rq->K = rq->src_part.IS;
			blklen = rq->src_part.IS * rq->T;
			test_log("SBN %u: short block %zu bytes\n", SBN, blklen);
		}
		if (SBN + 1 == rq->Z && rq->Kt * rq->T > rq->F) {
			/* last block needs padding */
			size_t padbyt = rq->Kt * rq->T - rq->F;
			test_log("last block and it needs padding of %zu bytes\n",
					padbyt);
			padblk = malloc(blklen);
			memcpy(padblk, srcblk, blklen - padbyt);
			memset(padblk + blklen - padbyt, 0, padbyt);
			srcblk = padblk;
		}

		rq_block(rq);
		rq_dump(rq, stderr);

		test_log("SBN %zu: K' (%u) * T (%u) = %zu\n", SBN, rq->KP, rq->T, rq->KP * rq->T);

		test_log("generating matrix A\n");
		rq_generate_matrix_A(rq, &A, rq->KP);

		test_log("matrix A done (%i x %i) %zu bytes\n", A.rows, A.cols, A.size);
		A_dup = matrix_dup(&A);
		test_assert(A.rows == A.cols, "Matrix A is square");
		test_assert(A.rows == rq->L, "Matrix A has dimension L");

		matrix_dump(&A, stderr);

		test_log("generating matrix D\n");
		D = rq_matrix_D(rq, srcblk);
		test_assert(D.rows = rq->L, "D has L rows");
		test_assert(D.cols = rq->T, "D has T cols");
		matrix_dump(&D, stderr);

		test_log("generating intermediate symbols\n");
		C = rq_intermediate_symbols(&A, &D);
		matrix_dump(&C, stderr);

		/* verify A*C=D */
		test_log("verifying A*C=D ... \n");
		matrix_t D2 = {0};
		matrix_multiply_gf256(&A_dup, &C, &D2);
		test_assert(D.rows == D2.rows, "rows of D match D2");
		test_assert(D.cols == D2.cols, "cols of D match D2");
		test_assert(memcmp(D.base, D2.base, D.size) == 0, "verify A*C=D");

		/* encode block */
		uint8_t *encsym, *repsym;
		encsym = malloc(rq->T * rq->K); assert(encsym);
		repsym = malloc(rq->T * nrep); assert(repsym);
		rq_encode_block(rq, &C, encsym, 0, rq->K);

		/* generate some repair symbols */
		rq_encode_block(rq, &C, repsym, rq->KP, nrep);

		matrix_t R = {0};
		matrix_new(&R, nrep, rq->T, repsym);
		fprintf(stderr, "Repair symbols:");
		matrix_dump(&R, stderr);

		decoder_tests(rq, srcblk, encsym, repsym);

		free(encsym);
		free(repsym);

		matrix_free(&D2);
		matrix_free(&D);
		matrix_free(&C);

		matrix_free(&A_dup);
		matrix_free(&A);

		size_t off = MIN(blklen, len);
		srcblk += off;
		len -= off;
		fprintf(stderr, "len remaining = %zu\n", len);
	}

	rq_free(rq);
	free(padblk);
	free(srcobj);

	test_log("test done\n");

	return fails;
}
