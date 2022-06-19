/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>
#include <sys/param.h>

#define REPS 1
#define FMIN 42
//#define FMAX 1017
#define FMAX 42
#define TMIN 8
#define TMAX 8
#define OMIN 5
#define OMAX 5

#define TEST_DEBUG 1

static_assert(TMIN % RQ_AL == 0);
static_assert(TMAX % RQ_AL == 0);

static uint32_t overhead;

uint8_t *generate_source_object(size_t F)
{
	uint8_t *obj = malloc(F);
	assert(obj);
	randombytes_buf(obj, F);
	test_log("object of %u bytes generated\n", F);
	return obj;
}

static int decodeC(rq_t *rq, uint8_t *enc, uint8_t *C1)
{
	uint8_t *ptr = C1;
	uint8_t *C2 = NULL;

	C2 = rq_decode_C(rq, enc);
	test_assert(C2 != NULL, "rq_decode_C()");

	fprintf(stderr, "C (original intermediate symbols):\n");
	for (uint32_t i = 0; i < rq->L; i++) {
		rq_dump_symbol(rq, ptr, stderr);
		ptr += rq->T;
	}

	if (C2) test_assert(!memcmp(C1, C2, rq->L * rq->T), "intermediate symbols match");

	free(C2);

	return 0;
}

static int bruteforce(rq_t *rq, matrix_t *A, uint8_t *enc, uint32_t ESI[], uint32_t nesi)
{
	uint8_t *dec;
	int rc = 0;

	dec = calloc(rq->K, rq->T);

	/* prepare Matrix A */
	rq_decoder_rfc6330_phase0(rq, A, dec, enc, ESI, nesi);

#if TEST_DEBUG
	fprintf(stderr, "Matrix A (%i x %i)\n", A->rows, A->cols);
	matrix_dump(A, stderr);
#endif
	matrix_gauss_elim(A, rq->sched);
	fprintf(stderr, "Schedule:\n");
	matrix_schedule_dump(rq->sched, stderr);

	free(dec);

	return rc;
}

static uint8_t *encoder_generate_symbols(rq_t *rq, uint32_t ESI[], int nesi)
{
	uint8_t *enc;
	uint8_t sbn = 0; /* FIXME - hardcoded SBN */
	rq_sym_t sym = {0};

	assert(rq->Z == 1); /* FIXME - only one block supported by this test */
	enc = malloc(nesi * rq->T);

	/* generate random repair symbols */
	sym.sym = enc;
	for (int i = 0; i < nesi; i++) {
		rq_symbol_random(rq, &sym, sbn); /* repair symbols, random order */
		ESI[i] = sym.ESI;
		sym.sym += rq->T;
	}
	return enc;
}

static int encoder_sizetest(uint8_t *src, size_t F, uint16_t T)
{
	rq_t *rq;
	matrix_t C;
	int nesi;
	int rc;

	rq = rq_init(F, T);
	rq_dump(rq, stderr);

	if (F == 0 || T == 0) {
		test_assert(!rq, "zero F or T, rq_init() returns  NULL");
		return 0;
	}
	rc = rq_encode_data(rq, src, F);
	test_assert(rc == 0, "rq_encode_data");

	C = rq_matrix_C_by_SBN(rq, 0);

	nesi = rq->KP + overhead; /* NB: overhead is over and above K', not K */
	if (!rc) {
		uint32_t ESI[nesi];
		uint8_t *enc = encoder_generate_symbols(rq, ESI, nesi);
		matrix_t A;
		rq_t *rq = rq_init(F, T); assert(rq);
		rq->sched = malloc(sizeof(matrix_sched_t));
		memset(rq->sched, 0, sizeof(matrix_sched_t));
		if (!rq->sched) {
			rq_free(rq);
			return -1;
		}

		assert(enc);

		rc = bruteforce(rq, &A, enc, ESI, nesi);
		test_assert(rc == 0, "Brute Force");
		rc = decodeC(rq, enc, C.base);
		test_assert(rc == 0, "Decode C");
		matrix_free(&A);
		rq_free(rq);
		free(enc);
	}
	rq_free(rq);
	return rc;
}

int main(void)
{
	uint8_t *srcobj = NULL;
	int rc;

	loginit();
	test_name("Hybrid (Gaussian + Decoding Schedule) Decoder");
	for (overhead = OMIN; overhead <= OMAX; overhead++) {
		for (int T = TMIN; T <= TMAX; T *= RQ_AL) {
			for (size_t F = FMIN; F <= FMAX; F++) {
				int ok = 0; int tests = 0;
				for (int i = 0; i < REPS; i++) {
					srcobj = generate_source_object(F);
					rc = encoder_sizetest(srcobj, F, T);
					free(srcobj);
					if (!rc) {
						test_log("SUCCESS with F=%zu, T=%u\n", F, T);
						ok++;
					}
					tests++;
					//if (rc) goto test_done; /* stop on first failure */
				}
				test_assert(ok == tests, "F=%zu, T=%u %i/%i (overhead = %u) tests ok",
						F, T, ok, tests, overhead);
			}
		}
	}
//test_done:
	test_log("test done\n");

	return fails;
}
