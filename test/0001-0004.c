/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <matrix.h>

int main(void)
{
	matrix_t A = {0};
	matrix_t LU = {0};
	matrix_t IA = {0};
	int P[4] = {0};
	int Q[4] = {0};

	loginit();
	test_name("Matrix LU Solve");

	uint8_t v0[] = {
		0x00, 0x11, 0x04, 0x05,
		0x0c, 0x01, 0x23, 0xef,
		0x02, 0x01, 0x03, 0x02,
		0x05, 0x09, 0x06, 0x0e
	};
	matrix_new(&A, 4, 4, v0, 0);

	LU = matrix_dup(&A);
	matrix_LU_decompose(&LU, P, Q);

	fprintf(stderr, "A:");
	matrix_dump(&A, stderr);
	fprintf(stderr, "LU:");
	matrix_dump(&LU, stderr);

	fprintf(stderr, "P[]:");
	for (int i = 0; i < (int)(sizeof P / sizeof P[0]); i++) {
		fprintf(stderr, " %i", P[i]);
	}
	fputc('\n', stderr);
	fprintf(stderr, "Q[]:");
	for (int i = 0; i < (int)(sizeof Q / sizeof Q[0]); i++) {
		fprintf(stderr, " %i", Q[i]);
	}
	fputc('\n', stderr);

	matrix_t B = matrix_dup(&A);
	matrix_zero(&B);
	matrix_solve_LU(&B, &A, &LU, P, Q);

	/* B is now the identity matrix */
	int cmp = 1;
	for (int i = 0; i < B.rows; i++) {
		for (int j = 0; j < B.cols; j++) {
			const uint8_t x = matrix_get(&B, i, j);
			const uint8_t y = (i == j) ? 1 : 0;
			if (x != y) cmp = 0;
		}
	}
	test_assert(cmp, "B is the identity matrix");

	fputc('\n', stderr);
	fprintf(stderr, "B:");
	matrix_dump(&B, stderr);

	matrix_free(&IA);
	matrix_free(&LU);
	matrix_free(&B);

	return fails;
}
