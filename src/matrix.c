/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <matrix.h>
#include <assert.h>
#include <stdlib.h>
#include <string.h>

matrix_t *matrix_new(matrix_t *mat, int rows, int cols, uint8_t *base)
{
	mat->rows = rows;
	mat->cols = cols;
	mat->trans = 0;
	mat->size = rows * cols;
	mat->base = (base) ? base : malloc(mat->size * sizeof(uint8_t));
	return mat;
}

matrix_t *matrix_zero(matrix_t *mat)
{
	memset(mat->base, 0, mat->size);
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
	size_t sz = m->size / m->cols / m->rows;
	assert(m->rows == m->cols);
	matrix_zero(m);
	for (int i = 0; i < m->cols * m->rows; i += m->cols + 1) {
		((char *)m->base)[i * sz] = 1;
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
	int r = (mat->trans) ? col : row;
	int c = (mat->trans) ? row : col;
	assert(r < mat->rows);
	assert(c < mat->cols);
	return mat->base[c + r * matrix_cols(mat)];
}

uint8_t matrix_set(matrix_t *mat, int row, int col, uint8_t val)
{
	int r = (mat->trans) ? col : row;
	int c = (mat->trans) ? row : col;
	assert(r < mat->rows);
	assert(c < mat->cols);
	mat->base[c + r * matrix_cols(mat)] = val;
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
		matrix_zero(p);
	}
	else {
		if (p->cols != x->rows || p->rows != y->rows)
			return NULL;
	}
	for (int i = 0; i < p->rows; i++) {
		for (int j = 0; j < p->cols; j++) {
			for (int k = 0; k < x->cols; k++) {
				const uint8_t a = matrix_get(x, i, k);
				const uint8_t b = matrix_get(y, k, j);
				matrix_inc_gf256(p, i, j, gf256_mul(a, b));
			}
		}
	}

	return p;
}

matrix_t *matrix_copy(matrix_t *dst, matrix_t *src)
{
	memcpy(dst->base, src->base, src->size);
	dst->rows = src->rows;
	dst->cols = src->cols;
	dst->trans = src->trans;
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
