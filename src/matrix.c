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
	mat->base = (base) ? base : malloc(mat->size);
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

void matrix_transpose(matrix_t *mat)
{
	mat->trans = !(mat->trans);
}

void matrix_free(matrix_t *mat)
{
	free(mat->base);
}
