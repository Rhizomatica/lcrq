/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <matrix.h>
#include <assert.h>
#include <gf256.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>

matrix_t *matrix_new(matrix_t *mat, const int rows, const int cols, uint8_t *base)
{
	mat->rows = rows;
	mat->cols = cols;
	mat->trans = 0;
	mat->stride = (size_t)cols * sizeof(uint8_t);
	mat->size = (size_t)rows * (size_t)cols * sizeof(uint8_t);
	mat->base = (base) ? base : malloc(mat->size * sizeof(uint8_t));
	return mat;
}

matrix_t matrix_submatrix(const matrix_t *A, const int off_rows, const int off_cols,
		const int rows, const int cols)
{
	const long int r = (A->trans) ? cols : rows;
	const long int c = (A->trans) ? rows : cols;
	const long int roff = (A->trans) ? off_cols : off_rows;
	const long int coff = (A->trans) ? off_rows : off_cols;
	assert(rows <= matrix_rows(A));
	assert(cols <= matrix_cols(A));
	matrix_t sub = {0};
	matrix_new(&sub, r, c, A->base + roff * A->stride + coff);
	sub.trans = A->trans;
	sub.stride = A->stride;
	return sub;
}

uint8_t *matrix_zero_row(matrix_t *m, int row)
{
	return memset(matrix_ptr_row(m, row), 0, m->cols);
}

/* A submatrix needs to be zeroed row by row.
 * A full matrix has a size which can be passed to memset */
void matrix_zero(matrix_t *m)
{
	if (!m->size) memset(m->base, 0, m->size);
	else for (int i = 0; i < m->rows; i++) matrix_zero_row(m, i);
}

matrix_t *matrix_identity(matrix_t *m)
{
	assert(m->rows == m->cols);
	matrix_zero(m);
	for (int i = 0; i < m->rows; i++) {
		matrix_set(m, i, i, 1);
	}
	return m;
}

/* return nonzero if m is an identity matrix */
int matrix_is_identity(matrix_t *m)
{
	int id = 1;
	for (int i = 0; i < m->rows; i++) {
		for (int j = 0; j < m->cols; j++) {
			const uint8_t x = matrix_get(m, i, j);
			const uint8_t y = (i == j) ? 1 : 0;
			if (x != y) id = 0;
		}
	}
	return id;
}

void matrix_dump(matrix_t *mat, FILE *stream)
{
	fprintf(stream, "\n");
	for (int r = 0; r < matrix_rows(mat); r++) {
		for (int c = 0; c < matrix_cols(mat); c++) {
			fprintf(stream, " %02x", matrix_get(mat, r, c));
		}
		fprintf(stream, "\n");
	}
	fprintf(stream, "\n");
}

void matrix_inc_gf256(matrix_t *mat, const int row, const int col, const uint8_t val)
{
	const uint8_t sum = matrix_get(mat, row, col) ^ val;
	matrix_set(mat, row, col, sum);
}

/* GF(256) dot product of x and y returned in p. Allocate p->base if required */
matrix_t *matrix_multiply_gf256(const matrix_t *x, const matrix_t *y, matrix_t *p)
{
	const int xrows = matrix_rows(x);
	const int yrows = matrix_rows(y);
	const int xcols = matrix_cols(x);
	const int ycols = matrix_cols(y);
	assert(xcols == yrows);
	if (!p->base) {
		matrix_new(p, xrows, ycols, NULL);
	}
	else {
		if (p->cols != xrows || p->rows != yrows)
			return NULL;
	}
	matrix_zero(p);
	for (int i = 0; i < p->rows; i++) {
		for (int j = 0; j < p->cols; j++) {
			uint8_t v = 0;
			for (int k = 0; k < xcols; k++) {
				const uint8_t a = matrix_get(x, i, k);
				const uint8_t b = matrix_get(y, k, j);
				v ^= GF256MUL(a, b);
			}
			matrix_set(p, i, j, v);
		}
	}

	return p;
}

matrix_t *matrix_swap_rows(matrix_t *m, const int r1, const int r2)
{
	for (int i = 0; i < matrix_cols(m); i++) {
		const uint8_t tmp = matrix_get(m, r1, i);
		matrix_set(m, r1, i, matrix_get(m, r2, i));
		matrix_set(m, r2, i, tmp);
	}
	return m;
}

matrix_t *matrix_swap_cols(matrix_t *m, const int c1, const int c2)
{
	for (int i = 0; i < matrix_rows(m); i++) {
		const uint8_t tmp = matrix_get(m, i, c1);
		matrix_set(m, i, c1, matrix_get(m, i, c2));
		matrix_set(m, i, c2, tmp);
	}
	return m;
}

void matrix_row_add(matrix_t *dst, const int drow, const matrix_t *src, const int srow)
{
	assert(matrix_cols(dst) == matrix_cols(src));
	uint8_t *dptr = matrix_ptr_row(dst, drow);
	const int mcols = matrix_cols(dst);
	for (int col = 0; col < mcols; col++) {
		*dptr ^= matrix_get(src, srow, col);
		dptr++;
	}
}

matrix_t matrix_add(const matrix_t *x, const matrix_t *y)
{
	assert(matrix_rows(x) == matrix_rows(y));
	assert(matrix_cols(x) == matrix_cols(y));
	matrix_t p = matrix_dup(x);
	matrix_zero(&p);
	for (int i = 0; i < matrix_rows(x); i++) {
		matrix_row_add(&p, i, x, i);
		matrix_row_add(&p, i, y, i);
	}
	return p;
}

void matrix_row_add_val(matrix_t *m, const int row, const uint8_t val)
{
	for (int col = 0; col < m->cols; col++) {
		matrix_set(m, row, col, matrix_get(m, row, col) ^ val);
	}
}

void matrix_row_div(matrix_t *m, const int row, const uint8_t val)
{
	for (int col = 0; col < m->cols; col++) {
		uint8_t a = matrix_get(m, row, col);
		uint8_t b = val;
		uint8_t v = GF256DIV(a, b);
		matrix_set(m, row, col, v);
	}
}

void matrix_col_mul(matrix_t *m, const int col, const int off, const uint8_t v)
{
	for (int row = off; row < matrix_rows(m); row++) {
		matrix_set(m, row, col, GF256MUL(matrix_get(m, row, col), v));
	}
}

void matrix_col_copy(matrix_t *dst, const int dcol, const matrix_t *src, const int scol)
{
	for (int row = 0; row < src->rows; row++) {
		matrix_set(dst, row, dcol, matrix_get(src, row, scol));
	}
}

void matrix_row_copy(matrix_t *dst, const int drow, const matrix_t *src, const int srow)
{
	uint8_t *dptr, *sptr;
	dptr = dst->base + drow * dst->stride;
	sptr = src->base + srow * src->stride;
	memcpy(dptr, sptr, src->stride);
}

static inline int pivot(matrix_t *A, int j, int P[], int Q[])
{
	const int Arows = matrix_rows(A);
	const int Acols = matrix_cols(A);
	for (int col = j; col < Acols; col++) {
		for (int row = j; row < Arows; row++) {
			if (matrix_get_s(A, row, j)) {
				/* pivot found, move in place, update P+Q */
				if (row != j) {
					matrix_swap_rows(A, row, j);
					SWAP(P[row], P[j]);
				}
				if (col != j) {
					matrix_swap_cols(A, col, j);
					SWAP(Q[col], Q[j]);
				}
				return 1;
			}
		}
	}
	return 0;
}

int matrix_LU_decompose(matrix_t *A, int P[], int Q[])
{
	const int Arows = matrix_rows(A);
	const int Acols = matrix_cols(A);
	const int n = MIN(Arows, Acols);
	int i;

	/* initialize permutations matricies */
	for (i = 0; i < Arows; i++) P[i] = i;
	for (i = 0; i < Acols; i++) Q[i] = i;

	/* LU decomposition */
	for (i = 0; i < n; i++) {
		uint8_t *Aii = MADDR(A, i, i);
		if (!pivot(A, i, P, Q)) break;
		for (int j = i + 1; j < Arows; j++) {
			uint8_t *Aji = MADDR(A, j, i);
			*Aji = GF256DIV(*Aji, *Aii);
			if (*Aji) for (int k = i + 1; k < Acols; k++) {
				const uint8_t b = matrix_get_s(A, i, k);
				A->base[k + j * A->stride] ^= GF256MUL(*Aji, b);
			}
		}
	}
	return i;
}

void matrix_inverse_LU(matrix_t *IA, const matrix_t *LU, const int P[])
{
	int n = MIN(matrix_rows(LU), matrix_cols(LU));

	assert(n == matrix_cols(LU));

	if (!IA->base) matrix_new(IA, matrix_rows(LU), matrix_cols(LU), NULL);

	for (int j = 0; j < matrix_cols(IA); j++) {
		for (int i = 0; i < matrix_cols(IA); i++) {
			uint8_t v = (P[i] == j) ? 1 : 0;
			matrix_set(IA, i, j, v);
			for (int k = 0; k < i; k++) {
				const uint8_t a = matrix_get(LU, i, k);
				const uint8_t b = matrix_get(IA, k, j);
				const uint8_t ab = GF256MUL(a, b);
				matrix_inc_gf256(IA, i, j, ab);
			}
		}

		for (int i = n - 1; i >= 0; i--) {
			for (int k = i + 1; k < matrix_cols(IA); k++) {
				const uint8_t a = matrix_get(LU, i, k);
				const uint8_t b = matrix_get(IA, k, j);
				const uint8_t ab = GF256MUL(a, b);
				matrix_inc_gf256(IA, i, j, ab);
			}
			const uint8_t v = gf256_div(matrix_get(IA, i, j), matrix_get(LU, i, i));
			matrix_set(IA, i, j, v);
		}
	}
}

void matrix_row_mul(matrix_t *m, const int row, const int off, const uint8_t val)
{
	for (int col = off; col < m->cols; col++) {
		matrix_set(m, row, col, GF256MUL(matrix_get(m, row, col), val));
	}
}

void matrix_row_mul_byrow(matrix_t *m, const int rdst, const int off, const int rsrc, const uint8_t factor)
{
	uint8_t *dptr = matrix_ptr_row(m, rdst) + off;
	uint8_t *sptr = matrix_ptr_row(m, rsrc) + off;
	for (int col = off; col < m->cols; col++) {
		if (*sptr && factor) {
			uint8_t f = GF256MUL(*sptr, factor);
			if (f) *dptr ^= f;
		}
		dptr++; sptr++;
	}
}

void matrix_solve_LU(matrix_t *X, const matrix_t *Y, const matrix_t *LU, const int P[], const int Q[])
{
	int n = MIN(matrix_rows(LU), matrix_cols(LU));

	if (!X->base) matrix_new(X, matrix_rows(LU), matrix_cols(Y), NULL);

	assert(matrix_cols(LU) == matrix_rows(X));
	assert(matrix_rows(LU) == matrix_rows(Y));
	assert(matrix_cols(X) == matrix_cols(Y));

	if (!X->base) matrix_new(X, matrix_rows(LU), matrix_cols(LU), NULL);

	for (int i = 0; i < n; i++) {
		matrix_row_copy(X, Q[i], Y, P[i]);
		for (int j = 0; j < i; j++) {
			const uint8_t LUij = matrix_get_s(LU, i, j);
			if (LUij) matrix_row_mul_byrow(X, Q[i], 0, Q[j], LUij);
		}
	}
	for (int i = n - 1; i >= 0; i--) {
		for (int j = i + 1; j < matrix_cols(LU); j++) {
			const uint8_t LUij = matrix_get_s(LU, i, j);
			if (LUij) matrix_row_mul_byrow(X, Q[i], 0, Q[j], LUij);
		}
		matrix_row_mul(X, Q[i], 0, GF256INV(matrix_get_s(LU, i, i)));
	}
}

matrix_t *matrix_inverse(matrix_t *A, matrix_t *I)
{
	matrix_new(I, A->rows, A->cols, NULL);
	matrix_identity(I);

	/* lower echelon form */
	for (int j = 0; j < A->cols; j++) {
		/* first, reduce the pivot row so jj = 1 */
		uint8_t jj = matrix_get(A, j, j);
		assert(jj); /* without row swaps, cannot proceed */
		if (jj != 1) {
			matrix_row_div(A, j, jj);
			matrix_row_div(I, j, jj);
		}
		for (int i = j + 1; i < A->rows; i++) {
			/* add pivot row (j) * factor to row i so that ij == 0 */
			jj = matrix_get(A, j, j);
			const uint8_t ij = matrix_get(A, i, j);
			const uint8_t f = gf256_div(ij, jj);
			matrix_row_mul_byrow(A, i, 0, j, f);
			matrix_row_mul_byrow(I, i, 0, j, f);
		}
	}
	/* finish upper triangle */
	for (int i = 0; i < A->cols; i++) {
		for (int j = i + 1; j < A->cols; j++) {
			const uint8_t ij = matrix_get(A, i, j);
			const uint8_t jj = matrix_get(A, j, j);
			const uint8_t f = gf256_div(ij, jj);
			matrix_row_mul_byrow(A, i, 0, j, f);
			matrix_row_mul_byrow(I, i, 0, j, f);
		}
	}

	return I;
}

matrix_t *matrix_copy(matrix_t *dst, const matrix_t *src)
{
	memcpy(dst->base, src->base, src->size);
	dst->rows = src->rows;
	dst->cols = src->cols;
	dst->trans = src->trans;
	dst->stride = src->stride;
	dst->size = src->size;
	return dst;
}

matrix_t matrix_dup(const matrix_t *src)
{
	matrix_t m;
	m.base = malloc(src->size);
	matrix_copy(&m, src);
	return m;
}

void matrix_transpose(matrix_t *mat)
{
	mat->trans = !(mat->trans);
}

void matrix_free(matrix_t *mat)
{
	free(mat->base);
	mat->base = NULL;
}
