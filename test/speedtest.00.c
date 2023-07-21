/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <arpa/inet.h>
#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <time.h>
#include <unistd.h>

#define DECODER_RFC 1
#define DECODER_DEFAULT DECODER_RFC

#define NANO 1000000000

#define DEFAULT_REPS 1
#define DEFAULT_F 42
#define DEFAULT_T 8
#define DEFAULT_O 1

static uint32_t overhead = RQ_OVERHEAD;
static int encoder_type = DECODER_DEFAULT;
static int decoder_type = DECODER_DEFAULT;

uint8_t *generate_source_object(size_t F)
{
	uint8_t *obj = malloc(F);
	assert(obj);
	test_randombytes(obj, F);
	return obj;
}

#if 0
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
#endif
static uint8_t *encoder_generate_symbols(rq_t *rq, uint32_t ESI[], int nesi)
{
	uint8_t *enc, *sym;
	rq_pid_t pid = 0;

	assert(rq->Z == 1); /* FIXME - only one block supported by this test */
	enc = malloc(nesi * rq->T);

	/* generate random repair symbols */
	sym = enc;
	for (int i = 0; i < nesi; i++) {
		rq_symbol(rq, &pid, sym, RQ_RAND);
		ESI[i] = rq_pid2esi(pid);
		sym += rq->T;
	}
	return enc;
}

uint8_t *encoder(rq_t *rq, uint8_t *src, uint32_t *ESI, int nesi)
{
	int rc = rq_encode(rq, src, rq->F);
	assert(rc == 0); (void)rc;
	return encoder_generate_symbols(rq, ESI, nesi);
}

static void dump_stats(const char *msg, size_t F, uint16_t T, size_t bytes, double s)
{
	fprintf(stderr, "%s %zu bytes in %0.4fs ", msg, bytes, s);
	double eBs = bytes / s;
	double ebps = eBs * 8;
	double eKbps = ebps / 1000;
	double eMbps = eKbps / 1000;
	double eGbps = eMbps / 1000;

	fprintf(stderr, "%0.1f Mbps, ", eMbps);
	fprintf(stderr, "%0.1f Gbps", eGbps);
	fprintf(stderr, " F=%zu, T=%u\n", F, T);
}

int main(int argc, char *argv[])
{
	rq_t *rq_enc = NULL, *rq_dec = NULL;
	double s_total_decoder = 0;
	double s_total_encoder = 0;
	struct timespec ts_enc_start = {0};
	struct timespec ts_enc_end = {0};
	struct timespec ts_dec_start = {0};
	struct timespec ts_dec_end = {0};
	uint8_t *srcobj, *enc;
	size_t bytes_total_decoder = 0;
	size_t bytes_total_encoder = 0;
	size_t F = DEFAULT_F;
	uint16_t T = DEFAULT_T;
	int reps = DEFAULT_REPS;
	uint32_t *ESI;
	int nesi;
	int fails = 0;

	if (argc > 1) F = atoll(argv[1]);
	if (argc > 2) T = atoll(argv[2]); /* TODO ensure multiple of Al */
	if (argc > 3) reps = atoll(argv[3]);
	if (argc > 4 && !strcmp(argv[4], "rfc")) {
		encoder_type = DECODER_RFC;
		decoder_type = DECODER_RFC;
	}

	for (int i = 0; i < reps; i++) {
		/* generate random source block for test */
		srcobj = generate_source_object(F);

		/* encoder test */
		clock_gettime(CLOCK_REALTIME, &ts_enc_start);
		rq_enc = rq_init(F, T);
		nesi = rq_enc->KP + overhead;
		ESI = calloc(nesi, sizeof(uint32_t));
		enc = encoder(rq_enc, srcobj, ESI, nesi);
		rq_free(rq_enc);
		clock_gettime(CLOCK_REALTIME, &ts_enc_end);

		/* encoder stats */
		uint64_t ensec = (ts_enc_end.tv_sec * NANO + ts_enc_end.tv_nsec);
		ensec -= (ts_enc_start.tv_sec * NANO + ts_enc_start.tv_nsec);
		double edsec = (double)ensec / NANO;
		bytes_total_encoder += F;
		s_total_encoder += edsec;

		/* decoder test */
		clock_gettime(CLOCK_REALTIME, &ts_dec_start);
		rq_dec = rq_init(F, T);
		uint8_t *dec = NULL;
		size_t decsz = rq_dec->K * rq_dec->T;
		dec = malloc(decsz);
		memset(dec, 0, decsz);
		int ok = rq_decode(rq_dec, dec, enc, ESI, nesi);
		rq_free(rq_dec);
		clock_gettime(CLOCK_REALTIME, &ts_dec_end);

		if (ok == 0) ok = memcmp(dec, srcobj, F);
		free(dec);

		/* decoder stats */
		uint64_t dnsec = (ts_dec_end.tv_sec * NANO + ts_dec_end.tv_nsec);
		dnsec -= (ts_dec_start.tv_sec * NANO + ts_dec_start.tv_nsec);
		double ddsec = (double)dnsec / NANO;
		s_total_decoder += ddsec;
		if (ok == 0) bytes_total_decoder += F;
		else fails++;

		/* clean up */
		free(ESI);
		free(enc);
		free(srcobj);
	}
	fputc('\n', stderr);
	dump_stats("encoder avg", F, T, bytes_total_encoder, s_total_encoder);
	dump_stats("decoder avg", F, T, bytes_total_decoder, s_total_decoder);
	if (fails) fprintf(stderr, "decoder FAILs = %i/%i (%0.4f%%), overhead = %i\n", fails, reps,
		(double)fails/(double)reps * 100, RQ_OVERHEAD);

	return 0;
}
