/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <arpa/inet.h>
#include <assert.h>
#include <lcrq_pvt.h>

static uint64_t F = 16;
static uint16_t T = 4;

static void dump_buffer(const uint8_t *sym, FILE *stream)
{
	for (uint64_t i = 0; i < F; i++) {
		fprintf(stream, " %02x", sym[i]);
	}
	fputc('\n', stream);
}

int main(void)
{
	rq_t *rq;
	uint8_t data[F];
	uint8_t copy[F];
	uint8_t *sym;
	uint16_t K;
	rq_pid_t pid = 0;
	int rc;

	loginit();
	test_name("Intermediate Symbol + Original Symbol Regeneration");

	/* generate random data */
	memset(data, 0, F);
	memset(copy, 0, F);
	test_randombytes(data, F);
	test_assert(memcmp(data, copy, F), "ensure buffers not match");

	/* calculate intermediate symbols */
	rq = rq_init(F, T);
	rc = rq_encode(rq, data, F);
	test_assert(rc == 0, "Encoding Successful");

	/* now re-create our original data from the intermediate symbols */
	K = rq_K(rq);
	sym = copy;
	test_log("K=%u\n", K);

	for (uint16_t esi = 0; esi < K; esi++) {
		pid = rq_pidsetesi(pid, esi);
		rq_symbol(rq, &pid, sym, 0);
		sym += T;
	}
	rq_free(rq);

	/* log results */
	fprintf(stderr, "orig: "); dump_buffer(data, stderr);
	fprintf(stderr, "copy: "); dump_buffer(copy, stderr);

	/* and compare */
	test_assert(!memcmp(data, copy, F), "regenerated data matches original");

	return test_status;
}
