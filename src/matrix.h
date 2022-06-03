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
	size_t   stride;
	size_t   size;
	uint8_t *base;
} matrix_t;

matrix_t *matrix_new(matrix_t *mat, const int rows, const int cols, uint8_t *base);
matrix_t matrix_submatrix(const matrix_t *A, const int off_rows, const int off_cols,
		const int rows, const int cols);
void matrix_free(matrix_t *mat);
matrix_t *matrix_zero(matrix_t *mat);
matrix_t *matrix_identity(matrix_t *mat);
void matrix_dump(matrix_t *mat, FILE *stream);
uint8_t matrix_get(const matrix_t *mat, const int row, const int col);
uint8_t matrix_set(matrix_t *mat, const int row, const int col, const uint8_t val);
int matrix_cols(const matrix_t *mat);
int matrix_rows(const matrix_t *mat);
void matrix_col_copy(matrix_t *dst, const int dcol, const matrix_t *src, const int scol);
void matrix_row_copy(matrix_t *dst, const int drow, const matrix_t *src, const int srow);
void matrix_row_add(matrix_t *dst, const int drow, const matrix_t *src, const int srow);
void matrix_row_add_val(matrix_t *m, const int row, const uint8_t val);

matrix_t matrix_add(const matrix_t *x, const matrix_t *y);

/* increment element by val using GF(256) addition */
uint8_t matrix_inc_gf256(matrix_t *mat, const int row, const int col, const uint8_t val);

/* GF(256) dot product of x and y returned in p. Allocate p->base if required */
matrix_t *matrix_multiply_gf256(const matrix_t *x, const matrix_t *y, matrix_t *p);

/* swap rows/cols in place */
matrix_t *matrix_swap_cols(matrix_t *m, const int c1, const int c2);
matrix_t *matrix_swap_rows(matrix_t *m, const int r1, const int r2);

/* I = (A^^-1) - set I to the inverse of A, allocating if required */
matrix_t *matrix_inverse(matrix_t *A, matrix_t *I);

/* mark matrix as "transposed" without modifying data */
void matrix_transpose(matrix_t *mat);

/* copy base data of src matrix to dst and set other vals */
matrix_t *matrix_copy(matrix_t *dst, const matrix_t *src);

/* duplicate matrix, allocating base */
matrix_t matrix_dup(const matrix_t *src);

#endif /* MATRIX_H */
