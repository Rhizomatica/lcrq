/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>

void rq_generate_HDPC(rq_t *rq, matrix_t *A);

char fileHDPC[] = "0000-0010.K%u.HDPC.tmp.XXXXXX";
char testHDPC[] = "hdpc.k%u.txt";

static void compare_files(FILE *f0, FILE *f1, uint16_t K)
{
	struct stat sb[2];
	char *map[2];
	int fd[2];

	fd[0] = fileno(f0);
	fd[1] = fileno(f1);

	for (int i = 0; i < 2; i++) {
		fstat(fd[i], &sb[i]);
		map[i] = mmap(NULL, sb[i].st_size, PROT_READ, MAP_PRIVATE, fd[i], 0);
		assert(map[i]);
		test_log("%i:\n%.*s\n", i, (int)sb[i].st_size, map[i]);
	}

	test_assert(sb[0].st_size == sb[1].st_size,
			"LDPC sizes match (K = %u) %zu == %zu",
			K, sb[0].st_size, sb[1].st_size);
	test_assert(!memcmp(map[0], map[1], sb[0].st_size), "LDPC codes match (K = %u)", K);

	for (int i = 0; i < 2; i++) {
		munmap(map[i], sb[i].st_size);
	}
}

static void dump_hdpc(const rq_t *rq, const matrix_t *A, FILE *stream)
{
	matrix_t H;
	H = matrix_submatrix(A, rq->S, 0, rq->H, rq->L);
	for (int r = 0; r < rq->H; r++) {
		for (int c = 0; c < rq->L; c++) {
			const uint8_t v = matrix_get(&H, r, c);
			fprintf(stream, " %02x", (int)v);
		}
		fputc('\n', stream);
	}
}

static void check_HDPC(rq_t *rq, uint16_t K)
{
	matrix_t A = {0};
	char ftmp[64];
	char cmpname[64];
	FILE *fout, *fcmp;

	sprintf(ftmp, fileHDPC, K);
	sprintf(cmpname, testHDPC, K);
	test_assert(mkstemp(ftmp) != -1, "mkstemp()");
	fout = fopen(ftmp, "w+");

	rq->K = K;
	rq_block(rq);
	rq_dump(rq, stderr);

	matrix_new(&A, rq->L, rq->L, NULL, 0);
	matrix_zero(&A);
	rq_generate_HDPC(rq, &A);
	dump_hdpc(rq, &A, fout);
	fflush(fout);
	matrix_free(&A);

	fcmp = fopen(cmpname, "r");
	test_log("diff %s %s\n", cmpname, ftmp);
	compare_files(fout, fcmp, K);

	fclose(fout);
	fclose(fcmp);
}

int main(void)
{
	rq_t *rq;

	loginit();
	test_name("check HDPC codes");

	rq = rq_init(1500, 1024);

	check_HDPC(rq, 1);
	check_HDPC(rq, 60);
	check_HDPC(rq, 4242);
	//check_HDPC(rq, 56403); // FIXME segfault, as does the TvRQ tool

	rq_free(rq);

	test_log("test done\n");

	return fails;
}
