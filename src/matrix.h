/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#ifndef MATRIX_H
#define MATRIX_H 1

#include <assert.h>
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

/* matrix_new() flags */
#define MATRIX_VEC 0x0001

/* schedule is stored as sequential records of matrix_op_*_t
 * the lower 4 bits of uchar type indicate the type of the following record
 * the upper 4 bits are the type of the previous record, so it is possible to
 * traverse the records in reverse */
enum {
	MATRIX_OP_NOOP = 0x0, /* end of records */
	MATRIX_OP_ROW  = 0x1,
	MATRIX_OP_COL  = 0x2,
	MATRIX_OP_ADD  = 0x3,
	MATRIX_OP_MUL  = 0x4
};
extern uint8_t reclen[5];

typedef struct matrix_op_add_s {
	uint8_t type;
	uint8_t beta;
	uint16_t dst;
	uint16_t src;
	uint16_t off;
} matrix_op_add_t;

typedef struct matrix_op_mul_s {
	uint8_t type;
	uint8_t beta;
	uint16_t dst;
} matrix_op_mul_t;

typedef struct matrix_op_swap_s {
	uint8_t type;
	uint16_t a;
	uint16_t b;
} matrix_op_swap_t;

typedef union {
	matrix_op_add_t add;
	matrix_op_mul_t mul;
	matrix_op_swap_t swp;
} matrix_op_t;

typedef struct matrix_sched_s {
	uint8_t *base;
	uint8_t *last;
	size_t len;
	size_t ops;
} matrix_sched_t;

typedef struct matrix_s {
	uint8_t *base;
	size_t   stride;
	size_t   size;
	int      cvec; /* cols aligned to vector multiple */
	int      rows;
	int      cols;
	int      roff;
	int      coff;
	int      trans;
} matrix_t;

matrix_t *matrix_new(matrix_t *mat, const int rows, const int cols, uint8_t *base, int flags);
matrix_t matrix_submatrix(const matrix_t *A, const int off_rows, const int off_cols,
		const int rows, const int cols);
void matrix_free(matrix_t *mat);

uint8_t *matrix_schedule_init(matrix_sched_t *sched);
uint8_t *matrix_schedule_resize(matrix_sched_t *sched);
void matrix_schedule_free(matrix_sched_t *sched);
void matrix_schedule_dump(matrix_sched_t *sched, FILE *stream);
void matrix_schedule_replay(matrix_t *m, matrix_sched_t *sched);
void matrix_schedule_reverse(matrix_t *m, matrix_sched_t *sched);

void matrix_sched_row(matrix_sched_t *sched, uint16_t a, uint16_t b);
void matrix_sched_col(matrix_sched_t *sched, uint16_t a, uint16_t b);
void matrix_sched_add(matrix_sched_t *sched, uint16_t dst, uint16_t off, uint16_t src, uint8_t beta);
void matrix_sched_mul(matrix_sched_t *sched, uint16_t dst, uint8_t beta);

/* Zero matrix row. Ignores transposition. Returns pointer to row */
uint8_t *matrix_zero_row(matrix_t *m, int row);

void matrix_zero(matrix_t *mat);
matrix_t *matrix_identity(matrix_t *mat);

/* return nonzero if m is an identity matrix */
int matrix_is_identity(matrix_t *m);

/* return nonzero if all the entries above the main diagonal are zero */
int matrix_is_lower(matrix_t *m);

/* return nonzero if all m elements are zero */
int matrix_is_zero(matrix_t *m);

void matrix_dump(matrix_t *mat, FILE *stream);
void matrix_col_copy(matrix_t *dst, const int dcol, const matrix_t *src, const int scol);
uint8_t *matrix_row_copy(matrix_t *dst, const int drow, const matrix_t *src, const int srow);
void matrix_row_add(matrix_t *dst, const int drow, const matrix_t *src, const int srow);
void matrix_row_add_val(matrix_t *m, const int row, const uint8_t val);
void matrix_col_mul(matrix_t *m, const int col, const int off, const uint8_t v);
void matrix_row_mul(matrix_t *m, const int row, const int off, const uint8_t val);

void matrix_row_mul_byrow(matrix_t *m, const int rdst, const int off, const int rsrc,
		const uint8_t factor);

matrix_t matrix_add(const matrix_t *x, const matrix_t *y);

/* increment element by val using GF(256) addition */
void matrix_inc_gf256(matrix_t *mat, const int row, const int col, const uint8_t val);

matrix_t *matrix_multiply_inplace(const matrix_t *x, matrix_t *y);

/* GF(256) dot product of x and y returned in p. Allocate p->base if required */
matrix_t *matrix_multiply_gf256(const matrix_t *x, const matrix_t *y, matrix_t *p);

/* swap rows/cols in place */
matrix_t *matrix_swap_cols(matrix_t *m, const int c1, const int c2);
matrix_t *matrix_swap_rows(matrix_t *m, const int r1, const int r2);

int matrix_pivot(matrix_t *A, int j, int P[], int Q[]);

int matrix_gauss_elim(matrix_t *A, matrix_sched_t *sched);

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
