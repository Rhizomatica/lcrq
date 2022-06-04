/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <matrix.h>

int main(void)
{
	matrix_t A = {0};
	matrix_t LU = {0};
	int Pr[4] = {0};

	loginit();
	test_name("Matrix LU Decompose");

	/* A . (A^^-1) == I_A
	 * generate the inverse of a matrix, then verify that the dot product of
	 * the two matrices is the identity matrix */
	uint8_t v0[] = {
		0x00, 0x11, 0x04, 0x05,
		0x0c, 0x01, 0x23, 0xef,
		0x02, 0x01, 0x03, 0x02,
		0x05, 0x09, 0x06, 0x0e
	};
	matrix_new(&A, 4, 4, v0);

	LU = matrix_dup(&A);
	matrix_LU_decompose(&LU, Pr);

	fprintf(stderr, "A:");
	matrix_dump(&A, stderr);
	fprintf(stderr, "LU:");
	matrix_dump(&LU, stderr);

	fprintf(stderr, "Pr[]:");
	for (int i = 0; i < (int)(sizeof Pr / sizeof Pr[0]); i++) {
		fprintf(stderr, " %i", Pr[i]);
	}
	fputc('\n', stderr);

#if 0
	A_orig = matrix_dup(&A);
	matrix_new(&I_A, 4, 4, NULL);
	matrix_identity(&I_A);

	matrix_inverse(&A, &I);

	test_log("inverse done, A, I:\n");

	matrix_dump(&A, stderr);
	matrix_dump(&I, stderr);

	/* the inverse of A has the same dimensions */
	test_assert(A.rows == I.rows, "I row count matches");
	test_assert(A.cols == I.cols, "I col count matches");
	test_assert(A.size == I.size, "I size matches");

	matrix_dump(&A_orig, stderr);
	matrix_dump(&I, stderr);

	/* multiply A by it's inverse I into P. */
	matrix_multiply_gf256(&A_orig, &I, &P);

	/* verify P == I_A */
	test_assert(P.rows == I_A.rows, "I_A row count matches");
	test_assert(P.cols == I_A.cols, "I_A col count matches");
	test_assert(P.size == I_A.size, "I_A size matches");
	test_assert(!memcmp(P.base, I_A.base, P.size), "data matches");

	matrix_dump(&P, stderr);

	mrtrix_free(&I);
	matrix_free(&I_A);
	matrix_free(&P);
	matrix_free(&A_orig);
#endif
	matrix_free(&LU);

	return fails;
}
