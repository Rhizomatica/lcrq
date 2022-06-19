/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>
#include <sys/param.h>

int main(void)
{
	loginit();

	test_name("Decoding Schedule replay");

	/* create a matrix with some data */
	matrix_t A = {0};
	matrix_new(&A, 10, 10, NULL);
	uint8_t v;
	for (int i = 0; i < A.rows; i++) {
		for (int j = 0; j < A.cols; j++) {
			randombytes_buf(&v, 1);
			matrix_set(&A, i, j, v);
		}
	}
	/* deliberately zero one of the main diagonal elements to force row swap */
	matrix_set(&A, 0, 0, 0);
	matrix_dump(&A, stderr);

	/* duplicate matrix */
	matrix_t B = matrix_dup(&A);

	/* create a blank schedule */
	matrix_sched_t sched = {0};
	matrix_schedule_init(&sched);

	/* reduce matrix, with tracking schedule */
	matrix_gauss_elim(&A, &sched);
	fprintf(stderr, "Matrix B:\n");
	matrix_dump(&A, stderr);

	/* replay schedule on duplicate matrix */
	matrix_schedule_replay(&B, &sched);

	fprintf(stderr, "Matrix B:\n");
	matrix_dump(&B, stderr);
	test_assert(!memcmp(A.base, B.base, A.size), "A matches B");

	matrix_schedule_free(&sched);
	matrix_free(&A);
	matrix_free(&B);

	return fails;
}
