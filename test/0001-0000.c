/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <matrix.h>

int main(void)
{
	uint8_t v;
	matrix_t A, B, C = {0};

	loginit();
	test_name("Matrix Multiply (GF256)");

	uint8_t v0[] = {
		1,
		1
	};
	uint8_t v1[] = {2};
	matrix_new(&A, 2, 1, v0, 0);
	matrix_new(&B, 1, 1, v1, 0);
	matrix_multiply_gf256(&A, &B, &C);
	v = matrix_get(&C, 0, 0);
	test_assert(C.base != NULL, "C base allocated");
	test_assert(C.rows == 2, "C.rows == 2");
	test_assert(C.cols == 1, "C.cols == 1");
	test_assert(v == 2, "C[0,0] == %u", v);
	v = matrix_get(&C, 1, 0);
	test_assert(v == 2, "C[1,0] == %u", v);
	matrix_dump(&C, stderr);
	matrix_free(&C);

	uint8_t v2[] = {
		4, 5, 6
	};
	uint8_t v3[] = {
		1,
		2,
		3
	};
	matrix_new(&A, 1, 3, v2, 0);
	matrix_new(&B, 3, 1, v3, 0);
	matrix_multiply_gf256(&A, &B, &C);
	test_assert(C.rows == 1, "C.rows == 1");
	test_assert(C.cols == 1, "C.cols == 1");
	v = matrix_get(&C, 0, 0);
	test_assert(v == 4, "C[0,0] == %u", v);
	matrix_dump(&C, stderr);
	matrix_free(&C);

	return test_status;
}
