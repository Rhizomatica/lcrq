/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <matrix.h>

int main(void)
{
	matrix_t A, A_orig, I = {0}, I_A, P = {0};
	loginit();
	test_name("Matrix Inverse");

	/* A . (A^^-1) == I_A
	 * generate the inverse of a matrix, then verify that the dot product of
	 * the two matrices is the identity matrix */

	uint8_t v0[] = {
		0, 17, 4,  5,
		12, 1,19, 14,
		2,  1, 3,  2,
		5,  9, 6, 11
	};
	matrix_new(&A, 4, 4, v0);
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

	matrix_free(&I);
	matrix_free(&I_A);
	matrix_free(&P);
	matrix_free(&A_orig);

	return fails;
}
