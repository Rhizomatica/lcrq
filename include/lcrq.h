/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#ifndef LCRQ_H
#define LCRQ_H 1

#include <matrix.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

/* rounding integer division */
#define CEIL(x, y)  (((x) + ((y) - 1)) / (y))
#define FLOOR(x, y) ((x) / (y))

typedef struct part_s {
	size_t IL;
	size_t IS;
	size_t JL;
	size_t JS;
} part_t;

typedef struct rq_tuple_s {
	uint16_t d;
	uint16_t a;
	uint16_t b;
	uint16_t d1;
	uint16_t a1;
	uint16_t b1;
} rq_tuple_t;

typedef struct rq_s rq_t;

part_t rq_partition(size_t I, uint16_t J);
size_t rq_rand(size_t y, uint8_t i, size_t m);
int rq_deg(rq_t *rq, int v);
rq_tuple_t rq_tuple(rq_t *rq, size_t X);

void rq_generate_matrix_A(rq_t *rq, matrix_t *A);
matrix_t rq_matrix_D(rq_t *rq, unsigned char *blk);
matrix_t rq_intermediate_symbols(matrix_t *A, matrix_t *D);
uint8_t *rq_encode(rq_t *rq, matrix_t *C, size_t isi);

void rq_block(rq_t *rq); /* calculate params based on K */
rq_t *rq_init(size_t F, uint16_t T);
void rq_free(rq_t *rq);
void rq_dump(rq_t *rq, FILE *stream);
void rq_dump_ldpc(rq_t *rq, matrix_t *A, FILE *stream);
void rq_dump_hdpc(rq_t *rq, matrix_t *A, FILE *stream);
void rq_dump_symbol(rq_t *rq, uint8_t *sym, FILE *stream);

#endif /* LCRQ_PVT_H */
