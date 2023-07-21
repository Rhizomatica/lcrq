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
	const int row_off = 0;
	const int col_off = 0;

	loginit();
	test_name("Matrix Submatrix");

	matrix_new(&A, 3, 12, NULL, 0);
	matrix_zero(&A);
	B = matrix_submatrix(&A, row_off, col_off, 3, 3);

	for (int i = 0; i < B.rows; i++) {
		for (int j = 0; j < B.cols; j++) {
			matrix_set(&B, i, j, 42);
		}
	}

	test_assert(B.rows == 3, "B.rows = %i", B.rows);
	test_assert(B.cols == 3, "B.cols = %i", B.cols);

	matrix_identity(&B);

	uint8_t a, b;
	for (int i = 0; i < A.rows; i++) {
		for (int j = 0; j < A.cols; j++) {
			b = (i < row_off || j < col_off) ? 0 : matrix_get(&B, i, j);
			a = matrix_get(&A, i, j);
			test_assert(a == b, "A[%i,%i]=%i", i, j, a);
		}
	}

	matrix_dump(&B, stderr);
	matrix_dump(&A, stderr);

	matrix_free(&A);

	return test_status;
}
