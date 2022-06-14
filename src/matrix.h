/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#ifndef MATRIX_H
#define MATRIX_H 1

#include <stdio.h>
#include <stddef.h>
#include <stdint.h>

#define SWAP(a, b) (((a) ^= (b)), ((b) ^= (a)), ((a) ^= (b)))
#define MCOL(m) (((m)->trans) ? (m)->rows : (m)->cols)
#define MROW(m) (((m)->trans) ? (m)->cols : (m)->rows)
#define MADDR(M, row, col) (M)->base + (row) * (M)->stride + (col)
#define MSET(M, row, col) (M)->base[(row) * (M)->stride + (col)]

#define matrix_get_s(M, row, col) (M)->base[(col) + (row) * (M)->stride]
#define matrix_get_t(M, row, col) (M)->base[(row) + (col) * (M)->stride]
#define matrix_get(M, row, col) \
	((M)->trans) \
	? (M)->base[(row) + (col) * (M)->stride] \
	: (M)->base[(col) + (row) * (M)->stride]

#define matrix_set_s(M, row, col, val) (M)->base[(col) + (row) * (M)->stride] = (val)
#define matrix_set_t(M, row, col, val) (M)->base[(row) + (col) * (M)->stride] = (val)
#define matrix_set(M, row, col, val) \
	if ((M)->trans) \
		(M)->base[(row) + (col) * (M)->stride] = val; \
	else \
		(M)->base[(col) + (row) * (M)->stride] = val

#define matrix_ptr_row(M, row) (M)->base + (row) * (M)->stride

#define matrix_cols MCOL
#define matrix_rows MROW

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

/* Zero matrix row. Ignores transposition. Returns pointer to row */
uint8_t *matrix_zero_row(matrix_t *m, int row);

void matrix_zero(matrix_t *mat);
matrix_t *matrix_identity(matrix_t *mat);

/* return nonzero if m is an identity matrix */
int matrix_is_identity(matrix_t *m);

/* return nonzero if all m elements are zero */
int matrix_is_zero(matrix_t *m);

void matrix_dump(matrix_t *mat, FILE *stream);
void matrix_col_copy(matrix_t *dst, const int dcol, const matrix_t *src, const int scol);
void matrix_row_copy(matrix_t *dst, const int drow, const matrix_t *src, const int srow);
void matrix_row_add(matrix_t *dst, const int drow, const matrix_t *src, const int srow);
void matrix_row_add_val(matrix_t *m, const int row, const uint8_t val);
void matrix_col_mul(matrix_t *m, const int col, const int off, const uint8_t v);

matrix_t matrix_add(const matrix_t *x, const matrix_t *y);

/* increment element by val using GF(256) addition */
void matrix_inc_gf256(matrix_t *mat, const int row, const int col, const uint8_t val);

/* GF(256) dot product of x and y returned in p. Allocate p->base if required */
matrix_t *matrix_multiply_gf256(const matrix_t *x, const matrix_t *y, matrix_t *p);

/* swap rows/cols in place */
matrix_t *matrix_swap_cols(matrix_t *m, const int c1, const int c2);
matrix_t *matrix_swap_rows(matrix_t *m, const int r1, const int r2);

/* peform LU decomposition on matrix A, storing combined LU factors in LU and
 * row and col permutations in P + Q. Return matrix rank */
int matrix_LU_decompose(matrix_t *A, int P[], int Q[]);

void matrix_inverse_LU(matrix_t *IA, const matrix_t *LU, const int P[]);

void matrix_solve_LU(matrix_t *X, const matrix_t *Y, const matrix_t *LU,
		const int P[], const int Q[]);

/* I = (A^^-1) - set I to the inverse of A, allocating if required */
matrix_t *matrix_inverse(matrix_t *A, matrix_t *I);

/* mark matrix as "transposed" without modifying data */
void matrix_transpose(matrix_t *mat);

/* copy base data of src matrix to dst and set other vals */
matrix_t *matrix_copy(matrix_t *dst, const matrix_t *src);

/* duplicate matrix, allocating base */
matrix_t matrix_dup(const matrix_t *src);

#endif /* MATRIX_H */
