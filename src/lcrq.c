/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <lcrq_pvt.h>
#include <matrix.h>
#include <assert.h>
#include <gf256.h>
#include <sodium.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <time.h>

static int isprime(const int n)
{
	if (n <= 1) return 0;
	if (n % 2 == 0) return 0;
	for (int i = 3; i < n; i += 2) {
		if (n % i == 0) return 0;
	}
	return 1;
}

/* return number of bits set in bitmap (Hamming Weight / popcount) */
static unsigned int hamm(unsigned char *map, size_t len)
{
	unsigned int c = 0;
	while (len--) for (char v = map[len]; v; c++) v &= v - 1;
	return c;
}

/* convert ESI to ISI (5.3.1) */
static inline uint32_t esi2isi(const rq_t *rq, uint32_t esi)
{
	assert(esi <= RQ_ESI_MAX);
	return (esi < rq->K) ? esi : esi + rq->KP - rq->K;
}

/* KL(n) is the maximum K' value in Table 2 in Section 5.6 such
	that K' <= WS/(Al*(ceil(T/(Al*n))))
	NB: this will always return K'_max = 56403 unless we have
	a small working memory (WS) */
static uint64_t KL(const uint64_t WS, const uint16_t Al, const uint16_t T, const uint16_t n)
{
	uint64_t v;
	for (int i = T2LEN - 1; i >= 0 ; i--) {
		v = WS/(Al*(CEIL(T,(Al*n))));
		assert(v);
		if (T2[i].k <= v) return T2[i].k;
	}
	return 0;
}

/* The degree generator Deg[v] is defined as follows, where v is a non-
   negative integer that is less than 2^^20 = 1048576.  Given v, find
   index d in Table 1 such that f[d-1] <= v < f[d], and set Deg[v] =
   min(d, W-2).  Recall that W is derived from K' as described in
   Section 5.3.3.3. */
int rq_deg(const rq_t *rq, const int v)
{
	int d;

	assert(v >= 0); assert(v < (1 << 20));

	for (d = 1; d <= DEGMAX; d++) {
		if (DEG[d-1] <= v && v < DEG[d]) break;
	}
	return MIN(d, rq->W-2);
}

size_t rq_rand(const size_t y, const uint8_t i, const size_t m)
{
	const uint8_t x0 = (y + i) % (1 << 8);
	const uint8_t x1 = ((y >> 8) + i) % (1 << 8);
	const uint8_t x2 = ((y >> 16) + i) % (1 << 8);
	const uint8_t x3 = ((y >> 24) + i) % (1 << 8);
	assert(m); /* must be positive */
	return (V0[x0] ^ V1[x1] ^ V2[x2] ^ V3[x3]) % m;
}

static matrix_t rq_matrix_C_by_SBN(const rq_t *rq, uint8_t SBN)
{
	matrix_t C = {0};
	assert(rq->C);
	matrix_new(&C, rq->L, rq->T, rq->C + SBN * rq->T * rq->L);
	return C;
}

/* Encoding Symbol Generator (5.3.5.3)
   The encoding symbol generator produces a single encoding symbol as
   output (referred to as result), according to the following algorithm:

	result = C[b]

	For j = 1, ..., d-1 do
		b = (b + a) % W
		result = result + C[b]

	While (b1 >= P) do b1 = (b1+a1) % P1

	result = result + C[W+b1]

	For j = 1, ..., d1-1 do
		b1 = (b1 + a1) % P1
		While (b1 >= P) do b1 = (b1+a1) % P1
		result = result + C[W+b1]

	Return result
*/
uint8_t *rq_encode_symbol(const rq_t *rq, const matrix_t *C, const uint32_t isi, uint8_t *sym)
{
	rq_tuple_t tup = rq_tuple(rq, isi);
	uint32_t b = tup.b;
	uint32_t b1 = tup.b1;
	matrix_t R;

	matrix_new(&R, 1, rq->T, sym);
	matrix_zero(&R);
	matrix_row_copy(&R, 0, C, b);
	for (uint32_t j = 1; j < tup.d; j++) {
		b = (b + tup.a) % rq->W;
		matrix_row_add(&R, 0, C, b);
	}
	while (b1 >= rq->P) b1 = (b1 + tup.a1) % rq->P1;
	matrix_row_add(&R, 0, C, rq->W + b1);
	for (uint32_t j = 1; j < tup.d1; j++) {
		b1 = (b1 + tup.a1) % rq->P1;
		while (b1 >= rq->P) b1 = (b1 + tup.a1) % rq->P1;
		matrix_row_add(&R, 0, C, rq->W + b1);
	}
	return R.base;
}

uint8_t *rq_symbol_generate(const rq_t *rq, rq_sym_t *sym, uint8_t sbn, uint32_t esi)
{
	uint32_t isi = esi2isi(rq, esi);
	matrix_t C = rq_matrix_C_by_SBN(rq, sbn);
	sym->SBN = sbn;
	sym->ESI = esi;
	return rq_encode_symbol(rq, &C, isi, sym->sym);
}

uint8_t *rq_symbol_random(const rq_t *rq, rq_sym_t *sym, uint8_t sbn)
{
	/* NB: ESI is a 24-bit unsigned integer (3.2) */
	uint32_t esi = randombytes_uniform(RQ_ESI_MAX - rq->K) + rq->K;
	return rq_symbol_generate(rq, sym, sbn, esi);
}

/* TODO - pass in state (threads) */
uint8_t *rq_symbol_repair_next(const rq_t *rq, rq_sym_t *sym, uint8_t sbn)
{
	static uint8_t sbn_last;
	static uint32_t esi;
	if (sbn != sbn_last || !esi || esi > RQ_ESI_MAX) esi = rq->K;
	sbn_last = sbn;
	return rq_symbol_generate(rq, sym, sbn, esi++);
}

/* reverse the polarity */
uint8_t *rq_symbol_repair_prev(const rq_t *rq, rq_sym_t *sym, uint8_t sbn)
{
	static uint8_t sbn_last;
	static uint32_t esi = RQ_ESI_MAX;
	if (sbn != sbn_last || esi < rq->K) esi = RQ_ESI_MAX;
	sbn_last = sbn;
	return rq_symbol_generate(rq, sym, sbn, esi--);
}

uint8_t *rq_encode_block(const rq_t *rq, const matrix_t *C, uint8_t *blk, uint32_t from, uint32_t n)
{
	for (uint32_t isi = from; isi < (from + n); isi++) {
		rq_encode_symbol(rq, C, isi, blk);
		blk += rq->T;
	}
	return blk;
}

int rq_encode_data(rq_t *rq, uint8_t *data, size_t len)
{
	const size_t clen = rq->L * rq->T;
	size_t blklen, off;
	uint8_t *base;
	uint8_t *padblk = NULL;
	const uint8_t *srcblk = data;

	/* create storage for intermediate symbols (C) */
	rq->C = malloc(clen * rq->Z);
	memset(rq->C, 0, clen * rq->Z);
	base = rq->C;
	rq->obj = data;
	rq->sym = data;
	for (uint8_t SBN = 0; SBN < rq->Z; SBN++) {
		matrix_t A, D;
		if (SBN < rq->src_part.JL) {
			rq->K = rq->src_part.IL;
			blklen = rq->src_part.IL * rq->T;
		}
		else {
			rq->K = rq->src_part.IS;
			blklen = rq->src_part.IS * rq->T;
		}
		if (SBN + 1 == rq->Z && rq->Kt * rq->T > rq->F) {
			/* last block needs padding */
			size_t padbyt = rq->Kt * rq->T - rq->F;
			fprintf(stderr, "padding last block with %zu bytes\n", padbyt);
			padblk = malloc(blklen);
			memcpy(padblk, srcblk, blklen - padbyt);
			memset(padblk + blklen - padbyt, 0, padbyt);
			srcblk = padblk;
		}
		rq_generate_matrix_A(rq, &A, rq->KP);
		D = rq_matrix_D(rq, srcblk);
		rq_intermediate_symbols(&A, &D, base);
		matrix_free(&D);
		matrix_free(&A);
		free(padblk);
		base += clen;
		off = MIN(blklen, len);
		srcblk += off;
		len -= off;
	}
	return 0;
}

uint8_t *rq_symbol_next(rq_state_t *state, rq_sym_t *sym)
{
	rq_t *rq = state->rq;
	if (state->ESI < rq->K && (state->flags & RQ_SOURCE) != RQ_SOURCE)
		state->ESI = rq->K;
	if (state->ESI >= rq->K) {
		if ((state->flags & RQ_REPAIR) == RQ_REPAIR) {
			matrix_t C = rq_matrix_C_by_SBN(rq, state->SBN);
			if (!state->rep) state->rep = malloc(rq->T);
			rq_encode_symbol(rq, &C, esi2isi(rq, state->ESI), state->rep);
			sym->sym = state->rep;
		}
		else {
			if (state->SBN + 1 >= rq->Z) { /* last block reached */
				if ((state->flags & RQ_REPEAT) == RQ_REPEAT) {
					state->SBN = 0;
					state->ESI = 0;
				}
				else return NULL;
			}
			state->SBN++; state->ESI = 0; /* advance to next block */
		}
	}
	if (state->ESI < rq->K) {
		assert(state->ESI * rq->T <= rq->F);
		sym->sym = state->sym;
		state->sym += rq->T;
	}
	sym->SBN = state->SBN;
	sym->ESI = state->ESI++;
	return sym->sym;
}

void rq_state_free(rq_state_t *state)
{
	free(state->rep);
}

void rq_state_init(rq_t *rq, rq_state_t *state, int flags)
{
	memset(state, 0, sizeof(rq_state_t));
	state->rq = rq;
	state->flags = flags;
	state->sym = rq->obj;
}

rq_tuple_t rq_tuple(const rq_t *rq, const uint32_t X)
{
	rq_tuple_t tup = {0};
	const uint32_t A = (53591 + rq->J * 997) | 0x1;
	const uint32_t B = 10267 * (rq->J + 1);
	const uint32_t y = B + X * A; /* mod 2^^32 */
	const uint32_t v = rq_rand(y, 0, 1 << 20);

	tup.d = rq_deg(rq, v);
	tup.a = 1 + rq_rand(y, 1, rq->W-1);
	tup.b = rq_rand(y, 2, rq->W);
	tup.d1 = (tup.d < 4) ? 2 + rq_rand(X, 3, 2) : 2;
	tup.a1 = 1 + rq_rand(X, 4, rq->P1 - 1);
	tup.b1 = rq_rand(X, 5, rq->P1);

	/* a is a positive integer between 1 and W-1 inclusive */
	assert(tup.a >= 1 && tup.a < rq->W);
	/* b is a non-negative integer between 0 and W-1 inclusive */
	assert(tup.b < rq->W);
	/* d1 is a positive integer that has value either 2 or 3 */
	assert(tup.d1 == 2 || tup.d1 == 3);
	/* a1 is a positive integer between 1 and P1-1 inclusive */
	assert(tup.a1 >= 1 && tup.a1 < rq->P1);
	/* b1 is a non-negative integer between 0 and P1-1 inclusive */
	assert(tup.b1 < rq->P1);

	return tup;
}

/* The function Partition[I,J] derives parameters for partitioning a
   block of size I into J approximately equal-sized blocks.  More
   specifically, it partitions I into JL blocks of length IL and JS
   blocks of length IS.  The output of Partition[I, J] is the sequence
   (IL, IS, JL, JS), where IL = ceil(I/J), IS = floor(I/J), JL = I - IS
   * J, and JS = J - JL.  */
part_t rq_partition(const size_t I, const uint16_t J)
{
	part_t p = {0};
	assert(J > 0);
	p.IL = CEIL(I,J);
	p.IS = FLOOR(I,J);
	p.JL = I - p.IS * J;
	p.JS = J - p.JL;
	return p;
}

/* The first row of Matrix A consists of three sub-matrices:
 * G_LDPC1, the identity matrix I_S and G_LDPC2
 * See RFC 6330 (5.3.3.3) p23 */
void rq_generate_LDPC(const rq_t *rq, matrix_t *A)
{
	matrix_t L1, I_S;

	matrix_new(&L1, rq->S, rq->L, A->base);

	/* G_LDPC,1 (S x B) */
	for (int i = 0; i < rq->B; i++) {        // For i = 0, ..., B-1 do
		const int a = 1 + i / rq->S;     //   a = 1 + floor(i/S)
		int b = i % rq->S;               //   b = i % S
		matrix_set(&L1, b, i, 1);        //   D[b] = D[b] + C[i]
		b = (b + a) % rq->S;             //   b = (b + a) % S
		matrix_set(&L1, b, i, 1);        //   D[b] = D[b] + C[i]
		b = (b + a) % rq->S;             //   b = (b + a) % S
		matrix_set(&L1, b, i, 1);        //   D[b] = D[b] + C[i]
	}

	/* The identity matrix, I_S */
	I_S = matrix_submatrix(A, 0, rq->B, rq->S, rq->S);
	matrix_identity(&I_S);

	/* G_LDPC,2 (S x P) */
	for (int i = 0; i < rq->S; i++) {         // For i = 0, ..., S-1 do
		const int a = i % rq->P;          //   a = i % P
		const int b = (i + 1) % rq->P;    //   b = (i+1) % P
		matrix_set(&L1, i, rq->W + a, 1); //   D[i] = D[i] + C[W+a] + C[W+b]
		matrix_set(&L1, i, rq->W + b, 1);
	}
}

/* The second row of Matrix A has the HDPC codes followed by
 * the identity matrix I_H
 * See RFC 6330 (5.3.3.3) p25 */
void rq_generate_HDPC(const rq_t *rq, matrix_t *A)
{
	matrix_t H1, I_H;
	uint8_t val = 1;

	H1 = matrix_submatrix(A, rq->S, 0, rq->H, rq->L);
	I_H = matrix_submatrix(&H1, 0, rq->L - rq->H, rq->H, rq->H);

	for (int j = 0; j < rq->H; j++) {
		matrix_set(&H1, j, rq->KP + rq->S - 1, GF256EXP(j));
	}
	for (int j = rq->KP + rq->S - 2; j >= 0; j--) {
		for (int i = 0; i < rq->H; i++) {
			val = matrix_get(&H1, i, j + 1);
			val = GF256MUL(val, 2);
			matrix_set(&H1, i, j, val);
		}
		int a = rq_rand(j + 1, 6, rq->H);
		val = matrix_get(&H1, a, j);
		val ^= 1; /* GF add => XOR */
		matrix_set(&H1, a, j, val);
		a = (a + rq_rand(j + 1, 7, rq->H - 1) + 1) % rq->H;
		val = matrix_get(&H1, a, j);
		val ^= 1; /* GF add => XOR */
		matrix_set(&H1, a, j, val);
	}

	/* The identity matrix, I_H */
	matrix_identity(&I_H);
}

static void rq_generate_LT(const rq_t *rq, uint8_t lt[rq->L], uint32_t isi)
{
	rq_tuple_t tup = rq_tuple(rq, isi);

	lt[tup.b] = 1;
	for (uint32_t j = 1; j < tup.d; j++) {
		tup.b = (tup.b + tup.a) % rq->W;
		lt[tup.b] = 1;
	}
	while (tup.b1 >= rq->P) tup.b1 = (tup.b1 + tup.a1) % rq->P1;
	lt[rq->W + tup.b1] = 1;
	for (uint32_t j = 1; j < tup.d1; j++) {
		tup.b1 = (tup.b1 + tup.a1) % rq->P1;
		while (tup.b1 >= rq->P) tup.b1 = (tup.b1 + tup.a1) % rq->P1;
		lt[rq->W + tup.b1] = 1;
	}
}

static void rq_generate_GENC(const rq_t *rq, matrix_t *A)
{
	matrix_t LT;
	matrix_new(&LT, rq->KP, rq->L, A->base + (rq->S + rq->H) * rq->L);
	for (int isi = 0; isi < LT.rows; isi++) {
		rq_generate_LT(rq, matrix_ptr_row(&LT, isi), isi);
	}
}

void rq_generate_matrix_A(const rq_t *rq, matrix_t *A, uint32_t lt)
{
	matrix_new(A, rq->S + rq->H + lt, rq->L, NULL);
	matrix_zero(A);
	assert(rq->L == rq->KP + rq->S + rq->H); /* L = K'+S+H (5.3.3.3) */
	assert(rq->L == rq->W + rq->P);          /* L = W+P (5.3.3.3) */
	rq_generate_LDPC(rq, A);
	rq_generate_HDPC(rq, A);
	rq_generate_GENC(rq, A);
}

/* 5.3.3.4.2
 * D denote the column vector consisting of S+H zero symbols followed
	by the K' source symbols C'[0], C'[1], ..., C'[K'-1] */
matrix_t rq_matrix_D(const rq_t *rq, const unsigned char *blk)
{
	uint8_t *ptr;
	matrix_t D = {0};

	assert(rq->KP + rq->S + rq->H == rq->L);
	matrix_new(&D, rq->L, rq->T, NULL);
	matrix_zero(&D);
	/* first S + H symbols are zero */
	ptr = D.base + (rq->S + rq->H) * rq->T * sizeof(uint8_t);
	/* copy K symbols of size T into D */
	memcpy(ptr, blk, rq->K * rq->T);

	return D;
}

/* calculate intermediate symbols (C) such that:
 *   C = (A^^-1)*D
 * where:
 *   D denotes the column vector consisting of S+H zero symbols
 *   followed by the K' source symbols C'[0], C'[1], ..., C'[K'-1]
 *
 *   if base is not NULL, symbols will be written here
 */
matrix_t rq_intermediate_symbols(matrix_t *A, const matrix_t *D, uint8_t *base)
{
	matrix_t LU = {0};
	matrix_t C = {0};
	int P[matrix_rows(A)];
	int Q[matrix_cols(A)];
	int rank;

	if (base) matrix_new(&C, matrix_rows(A), matrix_cols(D), base);
	LU = matrix_dup(A);
	rank = matrix_LU_decompose(&LU, P, Q);
	if (rank >= LU.cols) {
		LU.rows = rank;
		matrix_solve_LU(&C, D, &LU, P, Q);
	}
	matrix_free(&LU);

	return C;
}

void *rq_intermediate_symbols_alloc(const rq_t *rq)
{
	return calloc(rq->L, rq->L);
}

void rq_decoding_matrix_A(rq_t *rq, matrix_t *A, rq_blkmap_t *sym, rq_blkmap_t *rep)
{
	uint32_t lt = 0;

	/* create Matrix A (Directors extended cut) */

	matrix_new(A, rq->S + rq->H + rq->Nesi, rq->L, NULL);
	matrix_zero(A);
	rq_generate_LDPC(rq, A);
	rq_generate_HDPC(rq, A);

	/* append LT rows */

	/* first, any source symbols we have + padding rows */
	for (uint32_t isi = 0; isi < rq->KP; isi++) {
		if (isi < rq->K && !isset(sym->map, isi)) continue;
		const int row = rq->S + rq->H + lt++;
		rq_generate_LT(rq, matrix_ptr_row(A, row), isi);
	}
	/* repair symbols (K ... K + nrep) */
	for (int i = 0; i < rq->nrep; i++) {
		const uint32_t isi = esi2isi(rq, rep->ESI[i]);
		const uint32_t row = rq->S + rq->H + lt++;
		rq_generate_LT(rq, matrix_ptr_row(A, row), isi);
	}
}

static int rq_decoding_schedule(rq_t *rq, matrix_t *A, int P[], int Q[],
		rq_blkmap_t *sym, rq_blkmap_t *rep)
{
	int rank;

	rq_decoding_matrix_A(rq, A, sym, rep);

	/* LU decompose A */
	if ((rank = matrix_LU_decompose(A, P, Q)) < rq->L) {
		matrix_free(A);
		fprintf(stderr, "matrix_LU_decompose() FAIL rank = %i/%u\n", rank, rq->L);
		return -1;
	}
	A->rows = rank; /* discard extraneous rows */

	return 0;
}

static void rq_pack_LT_symbols(rq_t *rq, matrix_t *D, int P[], rq_blkmap_t *sym, rq_blkmap_t *rep)
{
	uint32_t gap = 0;
	uint32_t off = rq->S + rq->H;
	matrix_new(D, rq->L, rq->T, NULL);
	matrix_zero(D);
	for (int drow = 0, r = 0; drow < rq->KP; drow++) {
		const uint32_t srow = P[drow];
		if (drow < rq->nsrc + rq->nrep) {
			uint8_t *rptr;
			if (drow < rq->nsrc) {
				while (!isset(sym->map, srow + gap)) gap++;
				rptr = sym->p + rq->T * (srow + gap);
			}
			else if (drow < rq->nsrc + rq->KP - rq->K) {
				continue;
			}
			else {  /* repair symbol */
				assert(r < (int)rq->nrep);
				rptr = rep->p + rq->T * r++;
			}
			memcpy(matrix_ptr_row(D, drow + off), rptr, rq->T);
		}
	}
}

static int rq_decode_intermediate_symbols(rq_t *rq, rq_blkmap_t *sym, rq_blkmap_t *rep, matrix_t *C)
{
	const uint32_t M = rq->S + rq->H + rq->Nesi;
	matrix_t A = {0}, D = {0};
	int P[M], Q[rq->L];

	memset(P, 0, M);
	memset(Q, 0, rq->L);

	/* Build decoding schedule */
	if (rq_decoding_schedule(rq, &A, P, Q, sym, rep)) return -1;

	/* copy LT symbols into matrix */
	rq_pack_LT_symbols(rq, &D, P, sym, rep);

	/* solve to find intermediate symbols */
	matrix_new(C, rq->L, rq->T, NULL);
	matrix_zero(C);
	matrix_solve_LU(C, &D, &A, P, Q);

	matrix_free(&D);
	matrix_free(&A);

	return 0;
}

/* symbols will be written to sym directly as received. If the full set have
 * arrived, no decoding is required and they are already in place. If we are
 * missing any symbols, these symbols are added to the repair symbols and
 * decoding begins. We solve to find the intermediate symbols, then generate the
 * missing source symbols. */
int rq_decode_block(rq_t *rq, rq_blkmap_t *sym, rq_blkmap_t *rep)
{
	matrix_t C = {0};

	rq->nsrc = hamm(sym->map, sym->len);
	rq->nrep = (rep) ? hamm(rep->map, rep->len) : 0;
	rq->Nesi = rq->nsrc + rq->nrep + rq->KP - rq->K;

	/* if we have all of the first K symbols, no decoding required */
	if (rq->nsrc >= rq->K) return 0;

	/* check if we have enough symbols (including repair symbols) */
	/* NB: overhead is over and above K', not K */
	if (rq->Nesi < rq->KP) return -1;

	/* generate intermediate symbols */
	if (rq_decode_intermediate_symbols(rq, sym, rep, &C) == -1) return -1;

	/* generate missing source symbols */
	for (int esi = 0; esi < rq->K; esi++) {
		if (!isset(sym->map, esi)) {
			rq_encode_symbol(rq, &C, esi, sym->p + rq->T * esi);
		}
	}

	matrix_free(&C);

	return 0;
}

int rq_decode_block_f(rq_t *rq, uint8_t *dec, uint8_t *enc, uint32_t ESI[], uint32_t nesi)
{
	const size_t maplen = howmany(rq->K, CHAR_BIT);
	unsigned char symmap[maplen];
	unsigned char repmap[maplen];
	rq_blkmap_t sym = { .map = symmap, .len = maplen, .p = dec };
	rq_blkmap_t rep = { .map = repmap, .len = maplen, .p = enc, .ESI = ESI };

	memset(symmap, 0, maplen);
	memset(repmap, 0, maplen);

	for (uint32_t i = 0; i < nesi; i++) setbit(repmap, i);

	return rq_decode_block(rq, &sym, &rep);
}

void rq_decoder_rfc6330_phase0(rq_t *rq, matrix_t *A, uint8_t *dec, uint8_t *enc, uint32_t ESI[],
		uint32_t nesi)
{
	const size_t maplen = howmany(rq->K, CHAR_BIT);
	unsigned char symmap[maplen];
	unsigned char repmap[maplen];
	rq_blkmap_t sym = { .map = symmap, .len = maplen, .p = dec };
	rq_blkmap_t rep = { .map = repmap, .len = maplen, .p = enc, .ESI = ESI };

	memset(symmap, 0, maplen);
	memset(repmap, 0, maplen);

	for (uint32_t i = 0; i < nesi; i++) setbit(repmap, i);

	rq->nsrc = 0;
	rq->nrep = nesi;
	rq->Nesi = nesi + rq->KP - rq->K;

	rq_decoding_matrix_A(rq, A, &sym, &rep);
}

/*
	Phase 1 (5.4.2.2)
	+-----------+-----------------+---------+
	|           |                 |         |
	|     I     |    All Zeros    |         |
	|           |                 |         |
	+-----------+-----------------+    U    |
	|           |                 |         |
	|           |                 |         |
	| All Zeros |       V         |         |
	|           |                 |         |
	|           |                 |         |
	+-----------+-----------------+---------+
	Figure 6: Submatrices of A in the First Phase
*/

/* quick hack - HDPC rows are the only ones with values > 1 */
static int is_HDPC(matrix_t *A, int row)
{
	for (int j = 0; j < A->cols; j++) {
		if (matrix_get_s(A, row, j) > 1) return -1;
	}
	return 0;
}

int rq_phase1_choose_row(matrix_t *A, int i, int u, int *r, int odeg[],
		unsigned char comp[], int cmax, size_t mapsz)
{
	int rp = INT_MAX;
	int row = A->rows;
	int rdex[A->rows];

	memset(rdex, 0, sizeof rdex);

	/* Let r be the minimum integer such that at least one row of A has
	 * exactly r nonzeros in V */
	for (int x = i; x < A->rows; x++) {
		int r_row = 0;
		if (is_HDPC(A, x)) continue; /* skip HDPC rows */
		if (!hamm(matrix_ptr_row(A, x), A->stride)) continue;
		for (int y = i; y < A->cols - u; y++) {
			if (matrix_get_s(A, x, y)) r_row++;
			//if (r_row > rp) break; /* too high */
		}
		if (r_row > 0 && r_row < rp) {
			/* choose row with minimum original degree */
			//if (odeg[row] > r_row) row = x;
			row = x;
			rp = r_row;
			fprintf(stderr, "selecting row %i with r=%i\n", row, rp);
		}
		rdex[x] = r_row;
		fprintf(stderr, "row %i, r = %i\n", x, r_row);
	}

	/* FIXME - if a different row is chosen, update r */

	/* If r = 2 and there is a row with exactly 2 ones in V, then
	 * choose any row with exactly 2 ones in V that is part of a
	 * maximum size component in the graph described above that is
	 * defined by V. */
	if (rp == 2) {
		fprintf(stderr, "searching %i components...\n", cmax);
		unsigned char *cv;
		unsigned int component_sz = 0, sz;
		for (int v = 0; v < cmax; v++) {
			cv = comp + v * mapsz;
			sz = hamm(cv, mapsz); /* bits set => component size */
			if (sz == 0) break;   /* last component reached */
			if (sz > component_sz) {
				/* larger component found, choose row */
				component_sz = sz;
				for (int x = i; x < cmax; x++) {
					if (rdex[x] == 2 && isset(cv, x)) {
						row = x;
						fprintf(stderr, "component override\n");
						break;
					}
				}
			}
			else {
				fprintf(stderr, "component size %i not larger than %i\n",
						sz, component_sz);
			}
		}
	}
	*r = rdex[row];

	fprintf(stderr, "row %i was very much chosen with r=%i\n", row, *r);
	return (row == A->rows) ? -1 : row;
}

static void dump_components(unsigned char comp[], int cmax, size_t mapsz)
{
	for (int i = 0; i < cmax; i++) {
		fprintf(stderr, "%02i: ", i);
		if (hamm(comp + i * mapsz, mapsz))
		for (int j = 0; j < cmax; j++) {
			if (isset(comp + i * mapsz, j))
				fputc('1', stderr);
			else
				fputc('0', stderr);
		}
		fputc('\n', stderr);
	}
	fputc('\n', stderr);
}

/* track components in graph using a bitmap */
static void rq_graph_components(matrix_t *A, unsigned char comp[], int cmax,
		size_t mapsz, int i, int u)
{
	memset(comp, 0, mapsz * cmax);
	/* we're only interested in rows with exactly two nonzero elements */
	for (int x = i; x < A->rows; x++) {
		int r = 0;
		uint8_t a = 0, b = 0;
		if (is_HDPC(A, x)) continue; /* skip HDPC rows */
		if (!hamm(matrix_ptr_row(A, x), A->stride)) continue; /* skip zeros */
		for (int y = i; y < A->cols - u; y++) {
			if (matrix_get_s(A, x, y)) {
				if (r++ == 0) a = y;
				else b = y;
			}
			if (r > 2) break; /* too high */
		}
		if (r == 2) { /* row with r == 2, add to component bitmap */
			/* the two nonzero elements (a,b) are vertices. Setting the
			 * bits lets us find the largest component(s) by popcount */
			fprintf(stderr, "row %i is an edge connecting %i to %i, adding\n", x, a, b);
			int c[2] = { -1, -1};
			int v;
			unsigned char *cv, *c0, *c1;
			for (v = 0; v < cmax; v++) {
				cv = comp + v * mapsz;
				if (!hamm(cv, mapsz)) break; /* last component */
				if (isset(cv, a)) {
					fprintf(stderr, "setting a\n");
					c[0] = v;
					setbit(cv, b);
					if (c[1] > 0) break;
				}
				if (isset(cv, b)) {
					fprintf(stderr, "setting b\n");
					c[1] = v;
					setbit(cv, a);
					if (c[0] > 0) break;
				}
			}
			if (c[0] == -1 && c[1] == -1) {
				fprintf(stderr, "new component\n");
				/* new component */
				setbit(cv, a);
				setbit(cv, b);
				continue;
			}
			if (c[0] == c[1]) continue;
			if (c[0] > 0 && c[1] > 0) {
				/* merge two components into the first */
				size_t cl, cg;
				if (c[0] < c[1]) {
					cl = c[0];
					cg = c[1];
				}
				else {
					cl = c[1];
					cg = c[0];
				}
				fprintf(stderr, "merging component %zu into %zu\n", cg, cl);
				c0 = comp + cl * mapsz;
				c1 = comp + cg * mapsz;
				cv = comp + v * mapsz;
				for (int last = v; last < cmax; last++) {
					if (!hamm(cv, mapsz)) break; /* last component */
					cv = comp + v * mapsz;
				}
				dump_components(comp, cmax, mapsz); // FIXME
				for (size_t x = 0; x < mapsz; x++) {
					fprintf(stderr, "%i | %i == %i\n", c0[x], c1[x], c0[x] | c1[x]);
					c0[x] |= c1[x];
					c1[x] = 0;
					/* swap deleted component to end */
					if (cv != c1) c1[x] = cv[x]; //SWAP(c1[x], cv[x]);
				}
				fprintf(stderr, "merged:\n");
				dump_components(comp, cmax, mapsz); // FIXME
			}
		}
	}
}

int rq_decoder_rfc6330_phase1(rq_t *rq, matrix_t *X, matrix_t *A, int *i, int *u)
{
	int odeg[A->rows + 1];
	int row, r;
	int cmax = A->rows - rq->P;
	const size_t mapsz = howmany(cmax, CHAR_BIT);
	unsigned char *comp = malloc(cmax * mapsz);

	*i = 0;
	*u = rq->P;
	*X = matrix_dup(A);

	/* save original degree of each row */
	for (int i = 0; i < A->rows; i++) odeg[i] = matrix_row_degree(A, i);
	odeg[A->rows] = INT_MAX; /* last entry simplifies loop in row chooser */

	while ((*i) < A->rows && (*i) < A->cols && (*i) + (*u) < rq->L) {
		rq_graph_components(A, comp, cmax, mapsz, *i, *u);
		dump_components(comp, cmax, mapsz);
		/* all entries of V are zero => FAIL */
		if ((row = rq_phase1_choose_row(A, *i, *u, &r, odeg, comp, cmax, mapsz)) == -1) return -1;
		assert(row >= *i);

		//fprintf(stderr, "row %i chosen with r=%i (i = %u, u = %i)\n", row, r, *i, *u);

		/* the first row of A that intersects V is exchanged with the
		 * chosen row so that the chosen row is the first row that
		 * intersects V. */
		if (*i != row) {
			//fprintf(stderr, "swapping row %i with %i\n", *i, row);
			matrix_swap_rows(A, *i, row);
			matrix_swap_rows(X, *i, row);
			SWAP(odeg[*i], odeg[row]);
		}

		fprintf(stderr, "row swaps done:\n");
		matrix_dump(A, stderr);

		/* The columns of A among those that intersect V are
		 * reordered so that one of the r nonzeros in the chosen row
		 * appears in the first column of V and so that the remaining
		 * r-1 nonzeros appear in the last columns of V.  The same row
		 * and column operations are also performed on the matrix X. */
		int col = *i, j, rr = 0;
		int Vmax = A->cols - *u;
		for (j = *i + 1; j < Vmax; j++) {
			if (matrix_get_s(A, *i, j)) { /* nonzero found */
				if (matrix_get_s(A, *i, *i)) {
					/* swap to end */
					if (!rr) rr++;
					if (col == *i) col = Vmax - 1;
					while (matrix_get_s(A, *i, col)) col--;
					if (col <= j) break;
				}
				matrix_swap_cols(A, col, j);
				matrix_swap_cols(X, col, j);
				if (r == ++rr) break;
			}
		}

		fprintf(stderr, "col swaps done (i = %u, u = %i, r = %i)\n", *i, *u, r);

		matrix_t tmp = matrix_submatrix(A, 0, 0, *i + 1, Vmax);
		matrix_dump(&tmp, stderr);

		tmp = matrix_submatrix(A, 0, 0, *i + 1, *i + 1);
		matrix_dump(&tmp, stderr);

		tmp = matrix_submatrix(A, 0, *i + 1, *i + 1, Vmax - *i - 1);
		matrix_dump(&tmp, stderr);

		matrix_dump(A, stderr);

		/* Then, an appropriate multiple of the chosen row is added
		 * to all the other rows of A below the chosen row that have a
		 * nonzero entry in the first column of V.  Specifically, if a
		 * row below the chosen row has entry beta in the first column
		 * of V, and the chosen row has entry alpha in the first column
		 * of V, then beta/alpha multiplied by the chosen row is added
		 * to this row to leave a zero value in the first column of V.
		 * */
		const uint8_t alpha = matrix_get_s(A, *i, *i);
		const int ip = *i;
		for (int x = ip + 1; x < A->rows; x++) {
			const uint8_t beta = matrix_get_s(A, x, ip);
			if (beta > 0u) {
				const uint8_t f = GF256DIV(beta, alpha);
				matrix_row_mul_byrow(A, x, ip, ip, f);
			}
		}

		fprintf(stderr, "multiples added:\n");
		matrix_dump(A, stderr);

		/* i is incremented by 1 and u is incremented by r-1 */
		(*i)++;
		(*u) += r - 1;
		fprintf(stderr, "i=%i, u=%i\n", *i, *u);
	}
	fprintf(stderr, "i=%i, u=%i, L = %u\n", *i, *u, rq->L);

	/* trim X  - this is Phase 2 */
	//X->rows = *i;
	//X->cols = *i;

	free(comp);

	return 0;
}

#if 0
int rq_decoder_rfc6330(rq_t *rq, uint8_t *dec, uint8_t *enc, uint32_t ESI[], uint32_t nesi)
{
	return 0;
}
#endif

void rq_dump_hdpc(const rq_t *rq, const matrix_t *A, FILE *stream)
{
	matrix_t H;
	H = matrix_submatrix(A, rq->S, 0, rq->H, rq->L);
	for (int r = 0; r < rq->H; r++) {
		for (int c = 0; c < rq->L; c++) {
			const uint8_t v = matrix_get(&H, r, c);
			fprintf(stream, " %02x", (int)v);
		}
		fputc('\n', stream);
	}
}

void rq_dump_ldpc(const rq_t *rq, const matrix_t *A, FILE *stream)
{
	for (int r = 0; r < rq->S; r++) {
		for (int c = 0; c < rq->L; c++) {
			switch (matrix_get(A, r, c)) {
			case 0: fputc('0', stream); break;
			case 1: fputc('1', stream); break;
			default:
				fputc('-', stream); break;
			}
		}
		fputc('\n', stream);
	}
}

void rq_dump_symbol(const rq_t *rq, const uint8_t *sym, FILE *stream)
{
	//fprintf(stream, "symbol (%p)\n", (void *)sym);
	const int symwidth = MIN(rq->T, rq->F);
	for (int i = 0; i < symwidth; i++) {
		fprintf(stream, " %02x", sym[i]);
	}
	fputc('\n', stream);
}

void rq_free(rq_t *rq)
{
	if (rq) {
		free(rq->C);
		free(rq);
	}
}

/* dump all rq_t values to stream */
void rq_dump(const rq_t *rq, FILE *stream)
{
	fprintf(stream, "%s\t= %zu\n", "F", rq->F);
	fprintf(stream, "%s\t= %zu\n", "WS", rq->WS);
	fprintf(stream, "%s\t= %zu\n", "Kt", rq->Kt);
	fprintf(stream, "%s\t= %zu\n", "kl", rq->kl);
	fprintf(stream, "%s\t= %u\n", "Al", rq->Al);
	fprintf(stream, "%s\t= %u\n", "T", rq->T);
	fprintf(stream, "%s\t= %u\n", "SSS", rq->SSS);
	fprintf(stream, "%s\t= %u\n", "SS", rq->SS);
	fprintf(stream, "%s\t= %u\n", "Nmax", rq->Nmax);
	fprintf(stream, "%s\t= %u\n", "Z", rq->Z);
	fprintf(stream, "%s\t= %u\n", "N", rq->N);
	fprintf(stream, "%s\t= %u\n", "K", rq->K);
	fprintf(stream, "%s\t= %u\n", "K'", rq->KP);
	fprintf(stream, "%s\t= %u\n", "J", rq->J);
	fprintf(stream, "%s\t= %u\n", "H", rq->H);
	fprintf(stream, "%s\t= %u\n", "S", rq->S);
	fprintf(stream, "%s\t= %u\n", "W", rq->W);
	fprintf(stream, "%s\t= %u\n", "L", rq->L);
	fprintf(stream, "%s\t= %u\n", "P", rq->P);
	fprintf(stream, "%s\t= %u\n", "P1", rq->P1);
	fprintf(stream, "%s\t= %u\n", "U", rq->U);
	fprintf(stream, "%s\t= %u\n", "B", rq->B);
}

/* set per-block values */
void rq_block(rq_t *rq)
{
	for (int i = 0; i < T2LEN; i++) {
		if (T2[i].k >= rq->K) {
			rq->KP = T2[i].k;
			rq->H = T2[i].h;
			rq->S = T2[i].s;
			rq->W = T2[i].w;
			rq->J = T2[i].j;
			break;
		}
	}

	/* 5.3.3.3.  Pre-Coding Relationships */
	rq->L = rq->KP + rq->S + rq->H;
	rq->P = rq->L - rq->W;
	rq->U = rq->P - rq->H;
	rq->B = rq->W - rq->S;

	/* P1 denotes the smallest prime number greater than or equal to P */
	rq->P1 = rq->P;
	while (!isprime(rq->P1)) rq->P1++;

	/* N is the minimum n=1, ..., Nmax such that ceil(Kt/Z) <= KL(n) */
	for (rq->N = 1; rq->N <= rq->Nmax; rq->N++) {
		if (CEIL(rq->Kt,rq->Z) <= KL(rq->WS, rq->Al, rq->T, rq->N)) break;
	}

	/* verify primes */
	assert(isprime(rq->P1));
	assert(isprime(rq->S));
	assert(isprime(rq->W));
}

rq_t *rq_init(const size_t F, const uint16_t T)
{
	rq_t *rq;

	if (!F || !T) return NULL;

	rq = malloc(sizeof(rq_t));
	memset(rq, 0, sizeof(rq_t));

	rq->F = F;
	/* TODO what is an appropriate size for WS ? */
	rq->WS = 1073741824; /* 1GiB */
	rq->Al = RQ_AL;
	rq->T = T;

	/* T MUST be a multiple of Al */
	assert(T % rq->Al == 0);

	rq->SSS = T; /* sub-symbol size */
	rq->SS = rq->SSS / rq->Al;
	rq->Nmax = rq->T/(rq->SS*rq->Al);
	rq->Kt = CEIL(rq->F,rq->T);
	rq->kl = KL(rq->WS, rq->Al, rq->T, rq->Nmax);
	rq->Z = CEIL(rq->Kt,rq->kl);

	/* N is the minimum n=1, ..., Nmax such that ceil(Kt/Z) <= KL(n) */
	for (rq->N = 1; rq->N <= rq->Nmax; rq->N++) {
		if (CEIL(rq->Kt,rq->Z) <= KL(rq->WS, rq->Al, rq->T, rq->N)) break;
	}
	assert(rq->N == 1); /* no handling of sub-blocks */
	rq->src_part = rq_partition(rq->Kt, rq->Z);
	rq->sub_part = rq_partition(rq->T / rq->Al, rq->N);

	rq->K = rq->src_part.IL;
	rq_block(rq);

	return rq;
}
