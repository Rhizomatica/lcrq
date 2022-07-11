/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <arpa/inet.h>
#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sys/param.h>

#define REPS 1000
#define FMIN 42
#define FMAX 42
#define TMIN 8
#define TMAX 8
#define OMIN 1
#define OMAX 5

#ifdef static_assert
static_assert(TMIN % RQ_AL == 0, "TMIN not divisible by Al parameter");
static_assert(TMAX % RQ_AL == 0, "TMAX not divisible by Al parameter");
#endif

static uint32_t overhead;

uint8_t *generate_source_object(size_t F)
{
	uint8_t *obj = malloc(F);
	assert(obj);
	test_randombytes(obj, F);
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
	memset(dec, 0, rq->K * rq->T);
	if (!memcmp(dec, src, F)) return -1; /* ensure data doesn't match before decode */
	rc = rq_decode(rq, dec, enc, ESI, nesi);
	if (!rc) rc = memcmp(dec, src, F);
	free(dec);
	rq_free(rq);

	return rc;
}

static uint8_t *encoder_generate_symbols(rq_t *rq, uint32_t ESI[], int nesi)
{
	uint8_t *enc, *sym;
	rq_pid_t pid = 0;

	assert(rq->Z == 1); /* only one block supported by this test */
	enc = malloc(nesi * rq->T);

	/* generate random repair symbols */
	sym = enc;
	for (int i = 0; i < nesi; i++) {
		rq_symbol(rq, &pid, sym, RQ_RAND);
		ESI[i] = ntohl(pid << 8);
		sym += rq->T;
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
	rc = rq_encode(rq, src, F);
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
	test_name("5.3.3.4 Intermediate Symbol Generation + 5.3.4 Encoding");
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
