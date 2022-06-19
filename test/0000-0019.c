/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>
#include <sys/param.h>

#define REPS 100
#define FMIN 328183
#define FMAX 328183
#define TMIN 1024
#define TMAX 1024
#define OMIN 1
#define OMAX 1

static_assert(TMIN % RQ_AL == 0);
static_assert(TMAX % RQ_AL == 0);

// FIXME - T varies when == UINT16_MAX (or just large)

static uint32_t overhead;

uint8_t *generate_source_object(size_t F)
{
	uint8_t *obj = malloc(F);
	assert(obj);
	randombytes_buf(obj, F);
	test_log("object of %u bytes generated\n", F);
	return obj;
}

int decode_and_verify(uint8_t *src, uint8_t *enc, uint32_t ESI[], uint32_t nesi, size_t F, uint16_t T)
{
	rq_t *rq; /* separate rq context to ensure no stale data */
	uint8_t *dec;
	int rc = 0;

	rq = rq_init(F, T); assert(rq);
	dec = malloc(rq->K * rq->T);
	rc = rq_decode_block_f(rq, dec, enc, ESI, nesi);
	if (!rc) rc = memcmp(dec, src, F);
	free(dec);
	rq_free(rq);

	return rc;
}

uint8_t *encoder_generate_symbols(rq_t *rq, uint32_t ESI[], int nesi)
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

int encoder_sizetest(uint8_t *src, size_t F, uint16_t T)
{
	rq_t *rq;
	int nesi;
	int rc;

	rq = rq_init(F, T);

	if (F == 0 || T == 0) {
		test_assert(!rq, "zero F or T, rq_init() returns  NULL");
		return 0;
	}
	rc = rq_encode_data(rq, src, F);
	test_assert(rc == 0, "rq_encode_data");

	nesi = rq->KP + overhead; /* NB: overhead is over and above K', not K */
	if (!rc) {
		uint32_t ESI[nesi];
		uint8_t *enc = encoder_generate_symbols(rq, ESI, nesi);
		rc = decode_and_verify(src, enc, ESI, nesi, F, T);
		//test_assert(rc == 0, "decode and verify");
		free(enc);
	}

	/* something failed, dump rq object to logs */
	//if (rc) rq_dump(rq, stderr);
	rq_dump(rq, stderr);

	rq_free(rq);
	return rc;
}

int main(void)
{
	uint8_t *srcobj = NULL;
	int rc;

	loginit();
	return test_skip("Encoder Size Tests");
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
