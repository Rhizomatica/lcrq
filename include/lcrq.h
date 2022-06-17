/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#ifndef LCRQ_H
#define LCRQ_H 1

#include <matrix.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>

#define RQ_AL 4

/* allow for K + RQ_OVERHEAD repair symbols */
#define RQ_OVERHEAD 5

/* NB: ESI is a 24-bit unsigned integer (3.2) */
#define RQ_ESI_MAX 0xffffff

/* rounding integer division */
#define CEIL(x, y)  (((x) + ((y) - 1)) / (y))
#define FLOOR(x, y) ((x) / (y))

enum {
	RQ_SOURCE = 1,
	RQ_REPAIR = 2,
	RQ_REPEAT = 4,
	RQ_ASYNC  = 8
};

typedef struct rq_s rq_t;

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

typedef struct rq_sym_s {
	uint8_t *sym; /* pointer to symbol */
	uint32_t OIN; /* Object ID Number */
	uint32_t ESI; /* External Symbol ID */
	uint8_t  SBN; /* Source Block Number */
} rq_sym_t;

typedef struct rq_state_s {
	rq_t    *rq;
	uint8_t  SBN;
	uint32_t ESI;
	uint8_t * sym;
	uint8_t * rep;
	int      flags;
} rq_state_t;

typedef struct rq_blkmap_s {
	uint8_t *p;
	unsigned char *map;
	size_t len;
	uint32_t *ESI;
} rq_blkmap_t;

part_t rq_partition(const size_t I, const uint16_t J);
size_t rq_rand(size_t y, uint8_t i, size_t m);
int rq_deg(const rq_t *rq, const int v);
rq_tuple_t rq_tuple(const rq_t *rq, const uint32_t X);

void rq_generate_matrix_A(const rq_t *rq, matrix_t *A, uint32_t lt);
matrix_t rq_matrix_D(const rq_t *rq, const unsigned char *blk);
matrix_t rq_intermediate_symbols(matrix_t *A, const matrix_t *D, uint8_t *base);
uint8_t *rq_encode_symbol(const rq_t *rq, const matrix_t *C, const uint32_t isi, uint8_t *sym);

uint8_t *rq_symbol_generate(const rq_t *rq, rq_sym_t *sym, uint8_t sbn, uint32_t esi);

/* return a random repair symbol for block SBN
 * intermediate symbols must have been generated already by a call to
 * rq_encode_data() */
uint8_t *rq_symbol_random(const rq_t *rq, rq_sym_t *sym, uint8_t sbn);

uint8_t *rq_symbol_repair_next(const rq_t *rq, rq_sym_t *sym, uint8_t sbn);
uint8_t *rq_symbol_repair_prev(const rq_t *rq, rq_sym_t *sym, uint8_t sbn);

/* generate n symbols, starting at ISI from.  pass in preallocated buffer blk.
 * ISIs >= K are repair symbols */
uint8_t *rq_encode_block(const rq_t *rq, const matrix_t *C, uint8_t *blk, uint32_t from, uint32_t n);

/* encode an object of size F using a symbol size of T.
 * Return encoded data in enc, with size enclen */
//int rq_encode_object(const uint8_t *obj, const size_t F, const uint16_t T,
//		uint8_t *enc, size_t *enclen);

void rq_state_init(rq_t *rq, rq_state_t *state, int flags);
void rq_state_free(rq_state_t *state);
uint8_t *rq_symbol_next(rq_state_t *state, rq_sym_t *sym);

int rq_encode_data(rq_t *rq, uint8_t *data, const size_t len);

int rq_decode_block(rq_t *rq, rq_blkmap_t *sym, rq_blkmap_t *rep);

int rq_decode_block_f(rq_t *rq, uint8_t *dec, uint8_t *enc, uint32_t ESI[], uint32_t nesi);

void rq_decoding_matrix_A(rq_t *rq, matrix_t *A, rq_blkmap_t *sym, rq_blkmap_t *rep);
void rq_decoder_rfc6330_phase0(rq_t *rq, matrix_t *A, uint8_t *dec, uint8_t *enc, uint32_t ESI[],
		uint32_t nesi);
int rq_decoder_rfc6330_phase1(rq_t *rq, matrix_t *X, matrix_t *A, int *i, int *u);
int rq_decoder_rfc6330_phase2(rq_t *rq, matrix_t *A, matrix_t *X, int *i, int *u);
int rq_decoder_rfc6330_phase3(rq_t *rq, matrix_t *A, matrix_t *X, int *i, int *u);
int rq_decoder_rfc6330(rq_t *rq, uint8_t *dec, uint8_t *enc, uint32_t ESI[], uint32_t nesi);

void rq_block(rq_t *rq); /* calculate params based on K */

rq_t *rq_init(const size_t F, const uint16_t T);
void rq_free(rq_t *rq);
void rq_dump(const rq_t *rq, FILE *stream);
void rq_dump_ldpc(const rq_t *rq, const matrix_t *A, FILE *stream);
void rq_dump_hdpc(const rq_t *rq, const matrix_t *A, FILE *stream);
void rq_dump_symbol(const rq_t *rq, const uint8_t *sym, FILE *stream);

#endif /* LCRQ_PVT_H */
