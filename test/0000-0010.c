/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <string.h>
#include <sys/stat.h>
#include <sys/mman.h>
#include "valgrind.h"

void rq_generate_LDPC(rq_t *rq, matrix_t *A);

char fileLDPC[] = "0000-0010.K%u.LDPC.tmp.XXXXXX";
char testLDPC[] = "ldpc.k%u.txt";

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

static void dump_ldpc(const rq_t *rq, const matrix_t *A, FILE *stream)
{
	for (int r = 0; r < rq->S; r++) {
		for (int c = 0; c < rq->L; c++) {
			switch (matrix_get(A, r, c)) {
			case 0: fputc('0', stream); break;
			case 1: fputc('1', stream); break;
			default:
				fputc('-', stream); break;
			}
		}
		fputc('\n', stream);
	}
}

static void check_LDPC(rq_t *rq, uint16_t K)
{
	matrix_t A = {0};
	char ftmp[64];
	char cmpname[64];
	FILE *fout, *fcmp;

	sprintf(ftmp, fileLDPC, K);
	sprintf(cmpname, testLDPC, K);
	test_assert(mkstemp(ftmp) != -1, "mkstemp()");
	fout = fopen(ftmp, "w+");

	rq->K = K;
	rq_block(rq);
	rq_dump(rq, stderr);

	matrix_new(&A, rq->L, rq->L, NULL, 0);
	matrix_zero(&A);
	rq_generate_LDPC(rq, &A);
	dump_ldpc(rq, &A, fout);
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
	test_name("5.3.3.3 check LDPC codes");

	rq = rq_init(1500, 1024);

	check_LDPC(rq, 1);
	check_LDPC(rq, 60);
	check_LDPC(rq, 4242);
	/* testing with K'max under valgrind is slow */
	if (!RUNNING_ON_VALGRIND) check_LDPC(rq, 56403);

	rq_free(rq);

	test_log("test done\n");

	return fails;
}
