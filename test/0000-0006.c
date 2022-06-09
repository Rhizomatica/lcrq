/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>
#include <sys/param.h>

void rq_generate_HDPC(rq_t *rq, matrix_t *A);

static const size_t MAX_PAYLOAD = 4; /* MAX_PAYLOAD must be at least Al=4 bytes */
static const size_t MAX_SRCOBJ = MAX_PAYLOAD * 10 + 1;

/* generate source object of data of random size and content up to max bytes */
unsigned char *generate_source_object(size_t max, size_t *F)
{
	unsigned char *block;
	size_t sz = 0;
	//size_t sz = randombytes_uniform(max);
	if (!sz) sz = max;
	block = malloc(sz);
	assert(block);
	memset(block, 0, sz);
	//randombytes_buf(block, sz);
	for (unsigned char i = 0; i < MAX_SRCOBJ; i++) block[i] = i + 42;
	*F = sz;
	return block;
}

/*
The first set of pre-coding relations, called LDPC relations, is
described below and requires that at the end of this process the set
of symbols D[0] , ..., D[S-1] are all zero:

Initialize the symbols D[0] = C[B], ..., D[S-1] = C[B+S-1].

For i = 0, ..., B-1 do
	a = 1 + floor(i/S)
	b = i % S
	D[b] = D[b] + C[i]
	b = (b + a) % S
	D[b] = D[b] + C[i]
	b = (b + a) % S
	D[b] = D[b] + C[i]

For i = 0, ..., S-1 do
	a = i % P
	b = (i+1) % P
	D[i] = D[i] + C[W+a] + C[W+b]
*/
static void verify_LDPC_relations(rq_t *rq, matrix_t *C)
{
	matrix_t D = {0};
	matrix_new(&D, rq->S, rq->T, NULL);

	test_log("verifying LDPC coding relations\n");
	for (int i = 0; i < rq->S; i++) {
		matrix_row_copy(&D, i, C, i + rq->B);
	}
	matrix_dump(&D, stderr);

	uint8_t a, b;
	for (int i = 0; i < rq->B; i++) {
		a = 1 + FLOOR(i, rq->S);
		b = i % rq->S;
		matrix_row_add(&D, b, C, i);
		b = (b + a) % rq->S;
		matrix_row_add(&D, b, C, i);
		b = (b + a) % rq->S;
		matrix_row_add(&D, b, C, i);
	}
	for (int i = 0; i < rq->S; i++) {
		a = i % rq->P;
		b = (i + 1) % rq->P;
		matrix_row_add(&D, i, C, rq->W + a);
		matrix_row_add(&D, i, C, rq->W + b);
	}
	matrix_dump(&D, stderr);

	/* all entries in D MUST be zero */
	int checksum = 0;
	for (int i = 0; i < D.rows; i++) {
		for (int j = 0; j < D.cols; j++) {
			if (matrix_get(&D, i, j)) {
				checksum++;
				break;
			}
		}
	}
	test_assert(checksum == 0, "verify LDPC relations");
	matrix_free(&D);
}

/*
 * The second set of relations among the intermediate symbols C[0], ..., C[L-1]
 * are the HDPC relations and they are defined as follows:

Let:

alpha denote the octet represented by integer 2 as defined in Section 5.7.

MT denote an H x (K' + S) matrix of octets, where for j=0, ..., K'+S-2,
the entry MT[i,j] is the octet represented by the integer 1
if i= Rand[j+1,6,H] or i = (Rand[j+1,6,H] + Rand[j+1,7,H-1] + 1) % H,
and MT[i,j] is the zero element for all other values of i,
and for j=K'+S-1, MT[i,j] = alpha^^i for i=0, ..., H-1.

GAMMA denote a (K'+S) x (K'+S) matrix of octets, where
	GAMMA[i,j] =
		alpha ^^ (i-j) for i >= j,
		0 otherwise.

Then, the relationship between the first K'+S intermediate symbols
C[0], ..., C[K'+S-1] and the H HDPC symbols C[K'+S], ..., C[K'+S+H-1]
is given by:

Transpose[C[K'+S], ..., C[K'+S+H-1]] + MT * GAMMA *
Transpose[C[0], ..., C[K'+S-1]] = 0,

where '*' represents standard matrix multiplication utilizing the
octet multiplication to define the multiplication between a matrix of
octets and a matrix of symbols (in particular, the column vector of
symbols), and '+' denotes addition over octet vectors.
 */
static void verify_HDPC_relations(rq_t *rq, matrix_t *C)
{
	matrix_t MT, GAMMA, CT, CT1, CT2, P0, P1 = {0}, P2 = {0};

	/* build the MT matrix */
	matrix_new(&MT, rq->H, rq->KP + rq->S, NULL);
	matrix_zero(&MT);

	for (int j = 0; j < rq->KP +rq->S - 1; j++) {
		const uint32_t a = rq_rand(j + 1, 6, rq->H);
		const uint32_t b = rq_rand(j + 1, 7, rq->H - 1);
		const uint32_t c = (a + b + 1) % rq->H;
		matrix_set(&MT, a, j, 1);
		matrix_set(&MT, c, j, 1);
	}
	for (int i = 0; i < rq->H; i++) {
		const uint8_t v = gf256_exp(i);
		matrix_set(&MT, i, rq->KP + rq->S - 1, v);
	}

	fprintf(stderr, "MT:\n");
	matrix_dump(&MT, stderr);

	/* build GAMMA */
	matrix_new(&GAMMA, rq->KP + rq->S, rq->KP + rq->S, NULL);
	matrix_zero(&GAMMA);
	uint8_t v = 1;
	for (int i = 0; i < GAMMA.rows; i++) {
		for (int j = i; j < GAMMA.cols; j++) {
			matrix_set(&GAMMA, j, j - i, v);
		}
		v = gf256_mul(v, 2);
	}
	fprintf(stderr, "GAMMA:\n");
	matrix_dump(&GAMMA, stderr);

	matrix_multiply_gf256(&MT, &GAMMA, &P1);
	fprintf(stderr, "P1 (MT * GAMMA):\n");
	matrix_dump(&P1, stderr);

	/* ok, lets create the HDPC matrix another way */
	matrix_t A, HDPC;
	matrix_new(&A, rq->L, rq->L, NULL);
	rq_generate_HDPC(rq, &A);
	HDPC = matrix_submatrix(&A, rq->S, 0, rq->H, rq->KP + rq->S);
	fprintf(stderr, "HDPC:\n");
	matrix_dump(&HDPC, stderr);

	/* Tranpose C */
	CT = matrix_dup(C);
	matrix_transpose(&CT);
	fprintf(stderr, "Tranpose[C]:\n");
	matrix_dump(&CT, stderr);

	/* carve out submatricies from Tranpose[C] */
	CT1 = matrix_submatrix(&CT, 0, rq->KP + rq->S, rq->T, rq->H);
	CT2 = matrix_submatrix(&CT, 0, 0, rq->T, rq->KP + rq->S);

	fprintf(stderr, "CT1 = Transpose[C[K'+S], ..., C[K'+S+H-1]]:\n");
	matrix_dump(&CT1, stderr);
	fprintf(stderr, "CT2 = Transpose[C[0], ..., C[K'+S-1]]:\n");
	matrix_dump(&CT2, stderr);

/* RFC says:
	Transpose[C[K'+S], ..., C[K'+S+H-1]] + MT * GAMMA *
	Transpose[C[0], ..., C[K'+S-1]] = 0,

but, MT x GAMMA has dimensions (H x K'+S)
and Transpose[C[0], ..., C[K'+S-1]] has (T x K'+S)

For multiplication the columns of the first matrix MUST == rows of the second,

so we need to multiply the untransposed C[0], ..., C[K'+S-1] (K'+S x T)

(H x K'+S) x (K'+S x T) = (H x T)

which gives us a matrix P2(H x T)

and *then* transpose it to (T x K'+S) so we can add it to

Transpose[C[K'+S], ..., C[K'+S+H-1]] (T x H)

so why did we bother transposing either matrix?

we could have just added the two untransposed (H x T) matrices...?

*/
	matrix_transpose(&CT2); // we (re)transpose here ready to multiply
	matrix_multiply_gf256(&P1, &CT2, &P2);
	matrix_transpose(&P2);  // just so we can transpose back again to add
	fprintf(stderr, "P2 = (MT * GAMMA) * CT2:\n");
	matrix_dump(&P2, stderr);

	fprintf(stderr, "CT1:\n");
	matrix_dump(&CT1, stderr);

	fprintf(stderr, "P0 (CT1 + P2):\n");
	P0 = matrix_add(&CT1, &P2);
	matrix_dump(&P0, stderr);

	/* verify == 0 */
	int checksum = 0;
	for (int i = 0; i < matrix_rows(&P0); i++) {
		for (int j = 0; j < matrix_cols(&P0); j++) {
			if (matrix_get(&P0, i, j)) {
				checksum++;
				break;
			}
		}
	}
	test_assert(checksum == 0, "verify HDPC relations");

	matrix_free(&P0);
	matrix_free(&P1);
	matrix_free(&P2);
	matrix_free(&CT);
	matrix_free(&MT);
	matrix_free(&GAMMA);
	matrix_free(&A);
}

/*
Let G_ENC be the K' x L matrix such that

	G_ENC * Transpose[(C[0], ..., C[L-1])] =
	Transpose[(C'[0],C'[1], ...,C'[K'-1])],
*/
void verify_GENC(rq_t *rq, matrix_t *C, matrix_t *A)
{
	matrix_t G_ENC, CT, P = {0};

	G_ENC = matrix_submatrix(A, rq->S + rq->H, 0, rq->KP, rq->L);

	CT = matrix_dup(C);
	/* RFC says transpose, but the dimensions are correct without */
	//matrix_transpose(&CT);

	fprintf(stderr, "G_ENC:");
	matrix_dump(&G_ENC, stderr);
	fprintf(stderr, "C: (intermediate symbols)");
	matrix_dump(&CT, stderr);

	matrix_multiply_gf256(&G_ENC, &CT, &P);
	fprintf(stderr, "C': (source symbols)");
	matrix_dump(&P, stderr);

	/* TODO compare with source block C' */

	matrix_free(&CT);
	matrix_free(&P);
}

/* as per 5.3.2, the original source symbols C' can be generated
 * using the encoding process where 0 < ISI < K'
 * we will use this to test the encoder before generating repair
 * symbols */
void verify_encoder(rq_t *rq, matrix_t *C, uint8_t *src)
{
	uint8_t *sym;
	for (uint32_t isi = 0; isi < rq->K; isi++) {
		test_log("encoding ISI %zu\n", isi);
		sym = rq_encode_symbol(rq, C, isi, NULL);
		test_assert(!memcmp(sym, src, rq->T), "verify ISI %zu", isi);

		rq_dump_symbol(rq, src, stderr);
		rq_dump_symbol(rq, sym, stderr);

		free(sym);
		src += rq->T;
	}
}

int main(void)
{
	rq_t *rq;
	matrix_t A = {0};
	matrix_t D = {0};
	matrix_t A_dup = {0};
	matrix_t C;
	unsigned char *srcobj = NULL;
	unsigned char *padblk = NULL;
	unsigned char *srcblk;
	size_t F, len; /* F requires 40 bits */
	uint8_t SBN;
	size_t blklen;

	loginit();
	test_name("5.3.3.4 Intermediate Symbol Generation + 5.3.4 Encoding");

	/* generate random source block for test */
	srcobj = generate_source_object(MAX_SRCOBJ, &F);
	test_log("block of %u bytes generated\n", F);

	rq = rq_init(F, MAX_PAYLOAD);
	rq_dump(rq, stderr);

	srcblk = srcobj;
	len = F;

	test_log("ZL = %zu\n", rq->src_part.JL);
	test_log("ZS = %zu\n", rq->src_part.JS);

	/* encode blocks */
	for (SBN = 0; SBN < rq->Z; SBN++) {
		if (SBN < rq->src_part.JL) {
			/* long block */
			rq->K = rq->src_part.IL;
			blklen = rq->src_part.IL * rq->T;
			test_log("SBN %u: long block %zu bytes\n", SBN, blklen);
		}
		else {
			/* short block */
			rq->K = rq->src_part.IS;
			blklen = rq->src_part.IS * rq->T;
			test_log("SBN %u: short block %zu bytes\n", SBN, blklen);
		}
		if (SBN + 1 == rq->Z && rq->Kt * rq->T > rq->F) {
			/* last block needs padding */
			size_t padbyt = rq->Kt * rq->T - rq->F;
			test_log("last block and it needs padding of %zu bytes\n",
					padbyt);
			padblk = malloc(blklen);
			memcpy(padblk, srcblk, blklen - padbyt);
			memset(padblk + blklen - padbyt, 0, padbyt);
			srcblk = padblk;
		}

		rq_block(rq);
		rq_dump(rq, stderr);

		test_log("SBN %zu: K' (%u) * T (%u) = %zu\n", SBN, rq->KP, rq->T, rq->KP * rq->T);

		test_log("generating matrix A\n");
		rq_generate_matrix_A(rq, &A, rq->KP);

		//rq_dump_ldpc(rq, &A, stderr);
		//rq_dump_hdpc(rq, &A, stderr);

		test_log("matrix A done (%i x %i) %zu bytes\n", A.rows, A.cols, A.size);
		A_dup = matrix_dup(&A);
		test_assert(A.rows == A.cols, "Matrix A is square");
		test_assert(A.rows == rq->L, "Matrix A has dimension L");

		matrix_dump(&A, stderr);

		test_log("generating matrix D\n");
		D = rq_matrix_D(rq, srcblk);
		test_assert(D.rows = rq->L, "D has L rows");
		test_assert(D.cols = rq->T, "D has T cols");

		matrix_dump(&D, stderr);

		test_log("generating intermediate symbols\n");
		C = rq_intermediate_symbols(&A, &D, NULL);

		matrix_dump(&C, stderr);

		/* verify A*C=D */
		test_log("verifying A*C=D ... \n");
		matrix_t D2 = {0};
		matrix_multiply_gf256(&A_dup, &C, &D2);
		test_assert(D.rows == D2.rows, "rows of D match D2");
		test_assert(D.cols == D2.cols, "cols of D match D2");
		test_assert(memcmp(D.base, D2.base, D.size) == 0, "verify A*C=D");

		test_log("D:");
		matrix_dump(&D, stderr);
		test_log("D2:");
		matrix_dump(&D2, stderr);

		/* verify 5.3.3.3. Pre-Coding Relationships */
		verify_LDPC_relations(rq, &C);
		verify_HDPC_relations(rq, &C);
		verify_GENC(rq, &C, &A_dup);

		/* encoding (5.3.4) */
		verify_encoder(rq, &C, srcblk);

		/* TODO: repair symbols */
		//for (size_t isi = rq->KP; isi < rq->KP + 5; isi++) {
		//}

		matrix_free(&D2);
		matrix_free(&D);
		matrix_free(&C);

		matrix_free(&A_dup);
		matrix_free(&A);

		size_t off = MIN(blklen, len);
		srcblk += off;
		len -= off;
		fprintf(stderr, "len remaining = %zu\n", len);
	}

	rq_free(rq);
	free(padblk);
	free(srcobj);

	test_log("test done\n");

	return fails;
}
