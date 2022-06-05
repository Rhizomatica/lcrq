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
	test_name("Matrix LU Decompose + Inverse");

	/* A . (A^^-1) == IA
	 * generate the inverse of a matrix, then verify that the dot product of
	 * the two matrices is the identity matrix. We have a zero on the main
	 * diagonal here (0,0) so row swaps are required. */
	uint8_t v0[] = {
		0x00, 0x11, 0x04, 0x05,
		0x0c, 0x01, 0x23, 0xef,
		0x02, 0x01, 0x03, 0x02,
		0x05, 0x09, 0x06, 0x0e
	};
	matrix_new(&A, 4, 4, v0);

	LU = matrix_dup(&A);
	matrix_LU_decompose(&LU, P, Q);

	fprintf(stderr, "A:");
	matrix_dump(&A, stderr);
	fprintf(stderr, "LU:");
	matrix_dump(&LU, stderr);

	fprintf(stderr, "Pr[]:");
	for (int i = 0; i < (int)(sizeof P / sizeof P[0]); i++) {
		fprintf(stderr, " %i", P[i]);
	}
	fputc('\n', stderr);

	/* invert using LU */
	IA = matrix_dup(&A);
	matrix_inverse_LU(&IA, &LU, P);

	fprintf(stderr, "\nA^^-1:");
	matrix_dump(&IA, stderr);

	/* verify A x (A^^-1) = identity matrix */
	matrix_t R = {0};

	matrix_multiply_gf256(&A, &IA, &R);

	fprintf(stderr, "R = A x A^^-1 = Identity Matrix:");
	matrix_dump(&R, stderr);

	int cmp = 1;
	for (int i = 0; i < R.rows; i++) {
		for (int j = 0; j < R.cols; j++) {
			const uint8_t x = matrix_get(&R, i, j);
			const uint8_t y = (i == j) ? 1 : 0;
			if (x != y) {
				cmp = 0;
				fprintf(stderr, "(%i, %i): x=%u, y=%u\n", i, j, x, y);
				break;
			}
		}
	}
	test_assert(cmp, "result is identity matrix");

	matrix_free(&IA);
	matrix_free(&R);
	matrix_free(&LU);

	return fails;
}
