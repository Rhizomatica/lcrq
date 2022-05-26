/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#ifndef MATRIX_H
#define MATRIX_H 1

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

typedef struct matrix_s {
	int      rows;
	int      cols;
	int      trans;
	size_t   size;
	uint8_t *base;
} matrix_t;

matrix_t *matrix_new(matrix_t *mat, int rows, int cols, uint8_t *base);
void matrix_free(matrix_t *mat);
matrix_t *matrix_zero(matrix_t *mat);
matrix_t *matrix_identity(matrix_t *mat);
void matrix_dump(matrix_t *mat, FILE *stream);
uint8_t matrix_get(matrix_t *mat, int row, int col);
uint8_t matrix_set(matrix_t *mat, int row, int col, uint8_t val);

/* increment element by val using GF(256) addition */
uint8_t matrix_inc_gf256(matrix_t *mat, int row, int col, uint8_t val);

/* GF(256) dot product of x and y returned in p. Allocate p->base if required */
matrix_t *matrix_multiply_gf256(matrix_t *x, matrix_t *y, matrix_t *p);

/* mark matrix as "transposed" without modifying data */
void matrix_transpose(matrix_t *mat);

/* copy base data of src matrix to dst and set other vals */
matrix_t *matrix_copy(matrix_t *dst, matrix_t *src);

/* duplicate matrix, allocating base */
matrix_t matrix_dup(matrix_t *src);

#endif /* MATRIX_H */
