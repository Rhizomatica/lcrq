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

#endif /* MATRIX_H */
