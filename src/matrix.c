/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <matrix.h>
#include <assert.h>
#include <gf256.h>
#include <stdlib.h>
#include <string.h>

matrix_t *matrix_new(matrix_t *mat, int rows, int cols, uint8_t *base)
{
	mat->rows = rows;
	mat->cols = cols;
	mat->trans = 0;
	mat->stride = (size_t)cols * sizeof(uint8_t);
	mat->size = (size_t)rows * (size_t)cols * sizeof(uint8_t);
	mat->base = (base) ? base : malloc(mat->size * sizeof(uint8_t));
	return mat;
}

matrix_t matrix_submatrix(matrix_t *A, int off_rows, int off_cols, int rows, int cols)
{
	assert(rows <= A->rows);
	assert(cols <= A->cols);
	matrix_t sub = {0};
	matrix_new(&sub, rows, cols, A->base + off_rows * A->stride + off_cols);
	sub.stride = A->stride;
	return sub;
}

matrix_t *matrix_zero(matrix_t *mat)
{
	for (int i = 0; i < mat->rows; i++) {
		memset(mat->base + i * mat->stride, 0, mat->cols);
	}
	return mat;
}

int matrix_cols(matrix_t *mat)
{
	return (mat->trans) ? mat->rows : mat->cols;
}

int matrix_rows(matrix_t *mat)
{
	return (mat->trans) ? mat->cols : mat->rows;
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

void matrix_dump(matrix_t *mat, FILE *stream)
{
	uint8_t *i = mat->base;
	fprintf(stream, "\n");
	for (int r = 0; r < mat->rows; r++) {
		for (int c = 0; c < mat->cols; c++) {
			if (mat->trans) {
				fprintf(stream, "%4u", i[mat->cols * c + r]);
			}
			else {
				fprintf(stream, "%4u", *(i++));
			}
		}
		fprintf(stream, "\n");
	}
}

uint8_t matrix_get(matrix_t *mat, int row, int col)
{
	long int r = (mat->trans) ? col : row;
	long int c = (mat->trans) ? row : col;
	assert(r < mat->rows);
	assert(c < mat->cols);
	return mat->base[c + r * mat->stride];
}

uint8_t matrix_set(matrix_t *mat, int row, int col, uint8_t val)
{
	long int r = (mat->trans) ? col : row;
	long int c = (mat->trans) ? row : col;
	assert(r < mat->rows);
	assert(c < mat->cols);
	mat->base[c + r * mat->stride] = val;
	return val;
}

uint8_t matrix_inc_gf256(matrix_t *mat, int row, int col, uint8_t val)
{
	uint8_t sum = gf256_add(matrix_get(mat, row, col), val);
	return matrix_set(mat, row, col, sum);
}

/* GF(256) dot product of x and y returned in p. Allocate p->base if required */
matrix_t *matrix_multiply_gf256(matrix_t *x, matrix_t *y, matrix_t *p)
{
	assert(x->cols == y->rows);
	if (!p->base) {
		matrix_new(p, x->rows, y->cols, NULL);
	}
	else {
		if (p->cols != x->rows || p->rows != y->rows)
			return NULL;
	}
	matrix_zero(p);
	for (int i = 0; i < p->rows; i++) {
		for (int j = 0; j < p->cols; j++) {
			uint8_t v = 0;
			for (int k = 0; k < x->cols; k++) {
				const uint8_t a = matrix_get(x, i, k);
				const uint8_t b = matrix_get(y, k, j);
				v = gf256_add(v, gf256_mul(a, b));
			}
			matrix_set(p, i, j, v);
		}
	}

	return p;
}

matrix_t *matrix_swap_rows(matrix_t *m, int r1, int r2)
{
	for (int i = 0; i < matrix_cols(m); i++) {
		const uint8_t tmp = matrix_get(m, r1, i);
		matrix_set(m, r1, i, matrix_get(m, r2, i));
		matrix_set(m, r2, i, tmp);
	}
	return m;
}

matrix_t *matrix_swap_cols(matrix_t *m, int c1, int c2)
{
	for (int i = 0; i < matrix_rows(m); i++) {
		const uint8_t tmp = matrix_get(m, i, c1);
		matrix_set(m, i, c1, matrix_get(m, i, c2));
		matrix_set(m, i, c2, tmp);
	}
	return m;
}

void matrix_row_add(matrix_t *m, int row, uint8_t val)
{
	for (int col = 0; col < m->cols; col++) {
		matrix_set(m, row, col, gf256_add(matrix_get(m, row, col), val));
	}
}

void matrix_row_mul(matrix_t *m, int row, uint8_t val)
{
	for (int col = 0; col < m->cols; col++) {
		matrix_set(m, row, col, gf256_mul(matrix_get(m, row, col), val));
	}
}

void matrix_row_div(matrix_t *m, int row, uint8_t val)
{
	for (int col = 0; col < m->cols; col++) {
		uint8_t a = matrix_get(m, row, col);
		uint8_t b = val;
		uint8_t v = gf256_div(a, b);
		matrix_set(m, row, col, v);
	}
}

void matrix_row_mul_byrow(matrix_t *m, int rdst, int rsrc, uint8_t factor)
{
	for (int col = 0; col < m->cols; col++) {
		uint8_t dv = matrix_get(m, rdst, col);
		uint8_t sv = matrix_get(m, rsrc, col);
		uint8_t f = gf256_mul(sv, factor);
		if (f != 0) {
			matrix_set(m, rdst, col, gf256_add(dv, f));
		}
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
		if (jj != 1) {
			matrix_row_div(A, j, jj);
			matrix_row_div(I, j, jj);
		}
		for (int i = j + 1; i < A->rows; i++) {
			/* add pivot row (j) * factor to row i so that ij == 0 */
			jj = matrix_get(A, j, j);
			const uint8_t ij = matrix_get(A, i, j);
			const uint8_t f = gf256_div(ij, jj);
			matrix_row_mul_byrow(A, i, j, f);
			matrix_row_mul_byrow(I, i, j, f);
		}
	}
	/* finish upper triangle */
	for (int i = 0; i < A->cols; i++) {
		for (int j = i + 1; j < A->cols; j++) {
			const uint8_t ij = matrix_get(A, i, j);
			const uint8_t jj = matrix_get(A, j, j);
			const uint8_t f = gf256_div(ij, jj);
			matrix_row_mul_byrow(A, i, j, f);
			matrix_row_mul_byrow(I, i, j, f);
		}
	}

	return I;
}

matrix_t *matrix_copy(matrix_t *dst, matrix_t *src)
{
	memcpy(dst->base, src->base, src->size);
	dst->rows = src->rows;
	dst->cols = src->cols;
	dst->trans = src->trans;
	dst->stride = src->stride;
	dst->size = src->size;
	return dst;
}

matrix_t matrix_dup(matrix_t *src)
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
