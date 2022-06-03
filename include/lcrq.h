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
	uint32_t d;
	uint32_t a;
	uint32_t b;
	uint32_t d1;
	uint32_t a1;
	uint32_t b1;
} rq_tuple_t;

typedef struct rq_s rq_t;

part_t rq_partition(const size_t I, const uint16_t J);
size_t rq_rand(size_t y, uint8_t i, size_t m);
int rq_deg(const rq_t *rq, const int v);
rq_tuple_t rq_tuple(const rq_t *rq, const uint32_t X);

void rq_generate_matrix_A(const rq_t *rq, matrix_t *A);
matrix_t rq_matrix_D(const rq_t *rq, const unsigned char *blk);
matrix_t rq_intermediate_symbols(matrix_t *A, const matrix_t *D);
uint8_t *rq_encode(const rq_t *rq, const matrix_t *C, const uint32_t isi);

void rq_block(rq_t *rq); /* calculate params based on K */
rq_t *rq_init(const size_t F, const uint16_t T);
void rq_free(rq_t *rq);
void rq_dump(const rq_t *rq, FILE *stream);
void rq_dump_ldpc(const rq_t *rq, const matrix_t *A, FILE *stream);
void rq_dump_hdpc(const rq_t *rq, const matrix_t *A, FILE *stream);
void rq_dump_symbol(const rq_t *rq, const uint8_t *sym, FILE *stream);

#endif /* LCRQ_PVT_H */
