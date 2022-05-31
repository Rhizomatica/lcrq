/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <lcrq_pvt.h>
#include <matrix.h>
#include <assert.h>
#include <gf256.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <time.h>

int isprime(int n)
{
	if (n <= 1) return 0;
	if (n % 2 == 0) return 0;
	for (int i = 3; i < n; i++) {
		if (n % i == 0) return 0;
	}
	return 1;
}

/* KL(n) is the maximum K' value in Table 2 in Section 5.6 such
	that K' <= WS/(Al*(ceil(T/(Al*n))))
	NB: this will always return K'_max = 56403 unless we have
	a small working memory (WS) */
uint64_t KL(uint64_t WS, uint16_t Al, uint16_t T, uint16_t n)
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
int rq_deg(rq_t *rq, int v)
{
	assert(v >= 0);
	assert(v <= (1 << 20));

	int d;

	for (d = 1; d <= 30; d++) {
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
	return ((V0[x0] ^ V1[x1] ^ V2[x2] ^ V3[x3]) % m);
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
uint8_t *rq_encode(rq_t *rq, matrix_t *C, size_t isi)
{
	rq_tuple_t tup = rq_tuple(rq, isi);
	uint16_t b = tup.b;
	uint16_t b1 = tup.b1;
	matrix_t R;
	matrix_new(&R, 1, rq->T, NULL);
	matrix_zero(&R);

	matrix_row_copy(&R, 0, C, b);
	for (int j = 1; j < tup.d; j++) {
		b = (b ^ tup.a) % rq->W;
		matrix_row_add(&R, 0, matrix_get(C, b, 0));
	}
	while (b1 >= rq->P) b1 = (b1 ^ tup.a1) % rq->P1;
	matrix_row_add(&R, 0, matrix_get(C, rq->W + b1, 0));
	for (int j = 1; j < tup.d1; j++) {
		b1 = (b1 ^ tup.a1) % rq->P1;
		while (b1 >= rq->P) b1 = (b1 ^ tup.a1) % rq->P1;
		matrix_row_add(&R, 0, matrix_get(C, rq->W + b1, 0));
	}
	return R.base;
}


rq_tuple_t rq_tuple(rq_t *rq, size_t X)
{
	rq_tuple_t tup = {0};
	size_t A = 53591 + rq->J * 997;
	if (A % 2 == 0) A++;
	size_t B = 10267 * (rq->J + 1);
	size_t y = (B + X * A) & 0xffff;
	size_t v = rq_rand(y, 0, 1048576);
	size_t d = rq_deg(rq, v);
	tup.a = 1 + rq_rand(y, 1, rq->W-1);
	tup.b = rq_rand(y, 2, rq->W);
	tup.d1 = (d < 4) ? 2 + rq_rand(X, 3, 2) : 2;
	tup.a1 = 1 + rq_rand(X, 4, rq->P1 - 1);
	tup.b1 = rq_rand(X, 5, rq->P1);
	return tup;
}

/* The function Partition[I,J] derives parameters for partitioning a
   block of size I into J approximately equal-sized blocks.  More
   specifically, it partitions I into JL blocks of length IL and JS
   blocks of length IS.  The output of Partition[I, J] is the sequence
   (IL, IS, JL, JS), where IL = ceil(I/J), IS = floor(I/J), JL = I - IS
   * J, and JS = J - JL.  */
part_t rq_partition(size_t I, uint16_t J)
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
void rq_generate_LDPC(rq_t *rq, matrix_t *A)
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
void rq_generate_HDPC(rq_t *rq, matrix_t *A)
{
	matrix_t H1, I_H;
	uint8_t val = 1;

	H1 = matrix_submatrix(A, rq->S, 0, rq->H, rq->L);
	I_H = matrix_submatrix(&H1, 0, rq->L - rq->H, rq->H, rq->H);

	for (int j = 0; j < rq->H; j++) {
		matrix_set(&H1, j, rq->KP + rq->S - 1, val);
		val = gf256_mul(val, 2);
	}
	for (int j = rq->KP + rq->S - 2; j >= 0; j--) {
		for (int i = 0; i < rq->H; i++) {
			val = matrix_get(&H1, i, j + 1);
			val = gf256_mul(val, 2);
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

static void rq_generate_LT(rq_t *rq, matrix_t *A)
{
	matrix_t LT;
	matrix_new(&LT, rq->KP, rq->L, A->base + (rq->S + rq->H) * rq->L);
	for (int row = 0; row < LT.rows; row++) {
		uint32_t X = row - rq->S - rq->H;
		rq_tuple_t tup = rq_tuple(rq, X);
		matrix_set(&LT, row, tup.b, 1);
		for (int j = 1; j < tup.d; j++) {
			tup.b = (tup.b + tup.a) % rq->W;
			matrix_set(&LT, row, tup.b, 1);
		}
		while (tup.b1 >= rq->P)
			tup.b1 = (tup.b1 + tup.a1) % rq->P1;
		matrix_set(&LT, row, rq->W + tup.b1, 1);
		for (int j = 1; j < tup.d1; j++) {
			tup.b1 = (tup.b1 + tup.a1) % rq->P1;
			while (tup.b1 >= rq->P)
				tup.b1 = (tup.b1 + tup.a1) % rq->P1;
			matrix_set(&LT, row, rq->W + tup.b1, 1);
		}
	}
}

void rq_generate_matrix_A(rq_t *rq, matrix_t *A, uint8_t *src, size_t len)
{
	matrix_new(A, rq->L, rq->L, NULL);
	matrix_zero(A);
	assert(rq->L == rq->KP + rq->S + rq->H); /* L = K'+S+H (5.3.3.3) */
	assert(rq->L == rq->W + rq->P);          /* L = W+P (5.3.3.3) */
	rq_generate_LDPC(rq, A);
	rq_generate_HDPC(rq, A);
	rq_generate_LT(rq, A);
}

matrix_t rq_matrix_D(rq_t *rq, unsigned char *blk)
{
	uint8_t *ptr;
	matrix_t D = {0};

	matrix_new(&D, rq->L, 1, NULL);
	matrix_zero(&D);
	/* copy K' symbols of size T into D */
	ptr = D.base + (rq->S + rq->H);
	memcpy(ptr, blk, rq->KP);

	return D;
}

/* calculate intermediate symbols (C) such that:
 *   C = (A^^-1)*D
 * where:
 *   D denotes the column vector consisting of S+H zero symbols
 *   followed by the K' source symbols C'[0], C'[1], ..., C'[K'-1]
 */
matrix_t rq_intermediate_symbols(matrix_t *A, matrix_t *D)
{
	matrix_t A_inv = {0};
	matrix_t C = {0};
	fprintf(stderr, "%li: calculating A^^-1\n", clock());
	matrix_inverse(A, &A_inv);
	fprintf(stderr, "%li: multiplying (A^^-1).D => C\n", clock());
	matrix_multiply_gf256(&A_inv, D, &C);
	matrix_free(&A_inv);
	return C;
}

void *rq_intermediate_symbols_alloc(rq_t *rq)
{
	fprintf(stderr, "Matrix A: allocating %zu bytes\n", (size_t)rq->L * rq->L);
	return calloc(rq->L, rq->L);
}

void rq_dump_hdpc(rq_t *rq, matrix_t *A, FILE *stream)
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

void rq_dump_ldpc(rq_t *rq, matrix_t *A, FILE *stream)
{
	for (int r = 0; r < rq->S; r++) {
		for (int c = 0; c < rq->L; c++) {
			switch (matrix_get(A, r, c)) {
			case 0:         fputc('0', stream); break;
			case 1:         fputc('1', stream); break;
			default:
					fputc('-', stream); break;
			}
		}
		fputc('\n', stream);
	}
}

void rq_free(rq_t *rq)
{
	free(rq);
}

/* dump all rq_t values to stream */
void rq_dump(rq_t *rq, FILE *stream)
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

rq_t *rq_init(size_t F, uint16_t T)
{
	rq_t *rq = malloc(sizeof(rq_t));
	memset(rq, 0, sizeof(rq_t));

	rq->F = F;
	/* TODO what is an appropriate size for WS ? */
	rq->WS = 1073741824; /* 1GiB */
	rq->Al = 4;
	rq->T = T;
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
