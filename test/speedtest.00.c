/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>
#include <string.h>
#include <sys/param.h>
#include <time.h>
#include <unistd.h>

#define DECODER_GAUSS 1
#define DECODER_RFC 2
#define DECODER_HYBRID 3
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
	randombytes_buf(obj, F);
	return obj;
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

uint8_t *encoder_gauss(rq_t *rq, uint8_t *src, uint32_t *ESI, int nesi)
{
	int rc = rq_encode_data(rq, src, rq->F);
	assert(rc == 0);
	//fprintf(stderr, "gauss ");
	return encoder_generate_symbols(rq, ESI, nesi);
}

uint8_t *encoder_rfc(rq_t *rq, uint8_t *src, uint32_t *ESI, int nesi)
{
	int rc = rq_encode_data_rfc(rq, src, rq->F);
	assert(rc == 0);
	//fprintf(stderr, "rfc ");
	return encoder_generate_symbols(rq, ESI, nesi);
}

uint8_t *encoder(rq_t *rq, uint8_t *src, uint32_t *ESI, int nesi)
{
	uint8_t *(*f)(rq_t *, uint8_t *, uint32_t *, int);
	switch (encoder_type) {
		case DECODER_GAUSS:
			f = encoder_gauss;
			break;
		case DECODER_RFC:
			f = encoder_rfc;
			break;
		default:
			f = encoder_rfc;
	}
	return f(rq, src, ESI, nesi);
}

int decoder_gauss(uint8_t *enc, uint8_t *src, size_t F, uint16_t T, uint32_t ESI[], uint32_t nesi)
{
	rq_t *rq;
	uint8_t *dec;
	int rc = 0;

	//fprintf(stderr, "gauss ");

	rq = rq_init(F, T); assert(rq);
	dec = malloc(rq->K * rq->T);
	rc = rq_decode_block_f(rq, dec, enc, ESI, nesi);
	if (!rc) rc = memcmp(dec, src, F);
	free(dec);
	rq_free(rq);

	return rc;
}

int decoder_rfc(uint8_t *enc, uint8_t *src, size_t F, uint16_t T, uint32_t ESI[], uint32_t nesi)
{
	rq_t *rq;
	uint8_t *dec;
	size_t decsz;
	int rc = 0;

	//fprintf(stderr, "rfc ");

	rq = rq_init(F, T); assert(rq);
	decsz = rq->K * rq->T;
	dec = malloc(decsz);
	memset(dec, 0, decsz);
	rc = rq_decode_block_rfc(rq, dec, enc, ESI, nesi);
	if (!rc) rc = memcmp(dec, src, F);
	free(dec);
	rq_free(rq);

	return rc;
}

int decoder_hybrid(uint8_t *enc, uint8_t *src, size_t F, uint16_t T, uint32_t ESI[], uint32_t nesi)
{
	rq_t *rq;
	uint8_t *dec;
	size_t decsz;
	int rc = 0;

	//fprintf(stderr, "hybrid ");

	rq = rq_init(F, T); assert(rq);
	decsz = rq->K * rq->T;
	dec = malloc(decsz);
	memset(dec, 0, decsz);
	rc = rq_decode_block_hybrid(rq, dec, enc, ESI, nesi);
	if (!rc) rc = memcmp(dec, src, F);
	free(dec);
	rq_free(rq);

	return rc;
}

int decoder(uint8_t *enc, uint8_t *src, size_t F, uint16_t T, uint32_t ESI[], uint32_t nesi)
{
	int (*f)(uint8_t *, uint8_t *, size_t, uint16_t, uint32_t *, uint32_t);
	switch (decoder_type) {
		case DECODER_GAUSS:
			f = decoder_gauss;
			break;
		case DECODER_RFC:
			f = decoder_rfc;
			break;
		case DECODER_HYBRID:
			f = decoder_hybrid;
			break;
		default:
			f = decoder_rfc;
	}
	return f(enc, src, F, T, ESI, nesi);
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

	if (argc > 1) F = atoll(argv[1]);
	if (argc > 2) T = atoll(argv[2]); /* TODO ensure multiple of Al */
	if (argc > 3) reps = atoll(argv[3]);
	if (argc > 4 && !strcmp(argv[4], "gauss")) {
		encoder_type = DECODER_GAUSS;
		decoder_type = DECODER_GAUSS;
	}
	if (argc > 4 && !strcmp(argv[4], "hybrid")) decoder_type = DECODER_HYBRID;
	if (argc > 4 && !strcmp(argv[4], "rfc")) {
		encoder_type = DECODER_RFC;
		decoder_type = DECODER_RFC;
	}

	for (int i = 0; i < reps; i++) {

		/* generate random source block for test */
		srcobj = generate_source_object(F);

		/* encoder test */
		clock_gettime(CLOCK_REALTIME, &ts_enc_start);
		rq_t *rq = rq_init(F, T);
		nesi = rq->KP + overhead;
		ESI = calloc(nesi, sizeof(uint32_t));
		enc = encoder(rq, srcobj, ESI, nesi);
		rq_free(rq);
		clock_gettime(CLOCK_REALTIME, &ts_enc_end);

		/* encoder stats */
		uint64_t ensec = (ts_enc_end.tv_sec * NANO + ts_enc_end.tv_nsec);
		ensec -= (ts_enc_start.tv_sec * NANO + ts_enc_start.tv_nsec);
		double edsec = (double)ensec / NANO;

		//dump_stats("encoder", F, T, F, edsec);
		bytes_total_encoder += F;
		s_total_encoder += edsec;

		/* decoder test */
		clock_gettime(CLOCK_REALTIME, &ts_dec_start);
		int ok = decoder(enc, srcobj, F, T, ESI, nesi);
		clock_gettime(CLOCK_REALTIME, &ts_dec_end);

		/* decoder stats */
		if (ok == 0) {
			uint64_t dnsec = (ts_dec_end.tv_sec * NANO + ts_dec_end.tv_nsec);
			dnsec -= (ts_dec_start.tv_sec * NANO + ts_dec_start.tv_nsec);
			double ddsec = (double)dnsec / NANO;
			//dump_stats("decoder", F, T, F, ddsec);
			bytes_total_decoder += F;
			s_total_decoder += ddsec;
		}
		else {
			fprintf(stderr, "decoder FAIL ");
		}

		/* clean up */
		free(ESI);
		free(enc);
		free(srcobj);
	}
	fputc('\n', stderr);
	dump_stats("encoder avg", F, T, bytes_total_encoder, s_total_encoder);
	dump_stats("decoder avg", F, T, bytes_total_decoder, s_total_decoder);

	return 0;
}
