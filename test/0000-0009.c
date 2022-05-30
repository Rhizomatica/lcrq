/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <matrix.h>

/*
 * Create matrices A and B, where B
 * is a submatrix (view) of A
 *  +-----------------------+-------+
 *  |                       |       |
 *  |        A              |  B    |
 *  |                       |       |
 *  +-----------------------+-------+
 */

int main(void)
{
	matrix_t A = {0}, B = {0};

	loginit();
	test_name("Matrix Submatrix");

	matrix_new(&A, 3, 12, NULL);
	matrix_zero(&A);
	B = matrix_submatrix(&A, 0, 9, 3, 3);

	for (int i = 0; i < B.rows; i++) {
		for (int j = 0; j < B.cols; j++) {
			matrix_set(&B, i, j, 42);
		}
	}

	test_assert(B.rows == 3, "B.rows = %i", B.rows);
	test_assert(B.cols == 3, "B.cols = %i", B.cols);

	matrix_identity(&B);

	for (int i = 0; i < A.rows; i++) {
		for (int j = 0; j < A.cols; j++) {
			uint8_t v =
			((i == 0 && j == 9) ||
			(i == 1 && j == 10) ||
			(i == 2 && j == 11)) ? 1 : 0;
			test_assert(matrix_get(&A, i, j) == v,
					"A[%i,%i] == %u", i, j, v);
		}
	}

	matrix_dump(&B, stderr);
	matrix_dump(&A, stderr);

	matrix_free(&A);

	return fails;
}
