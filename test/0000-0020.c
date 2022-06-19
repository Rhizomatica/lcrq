/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <gf256.h>
#include <lcrq.h>
#include <lcrq_pvt.h>
#include <sodium.h>
#include <sys/param.h>

#define REPS 1
#define FMIN 42
//#define FMAX 1017
#define FMAX 42
#define TMIN 8
#define TMAX 8
#define OMIN 5
#define OMAX 5

#define TEST_DEBUG 1

static_assert(TMIN % RQ_AL == 0);
static_assert(TMAX % RQ_AL == 0);

static uint32_t overhead;

uint8_t *generate_source_object(size_t F)
{
	uint8_t *obj = malloc(F);
	assert(obj);
	randombytes_buf(obj, F);
	test_log("object of %u bytes generated\n", F);
	return obj;
}

static int decodeC(rq_t *rq, uint8_t *enc, uint8_t *C1)
{
	uint8_t *C2 = NULL;

	fprintf(stderr, "Schedule:\n");
	matrix_schedule_dump(rq->sched, stderr);

	C2 = rq_decode_C(rq, enc);
	test_assert(C2 != NULL, "rq_decode_C()");

	fprintf(stderr, "C (original intermediate symbols):\n");
	uint8_t *ptr = C1;
	for (uint32_t i = 0; i < rq->L; i++) {
		rq_dump_symbol(rq, ptr, stderr);
		ptr += rq->T;
	}

	if (C2) test_assert(!memcmp(C1, C2, rq->L * rq->T), "intermediate symbols match");

	free(C2);

	return 0;
}

static int phase_3(rq_t *rq, matrix_t *A, matrix_t *X, int i, int u)
{
	int rc = 0;

	matrix_t U_upper = matrix_submatrix(A, 0, i, i, u);

	/* TODO the matrix X is multiplied with the submatrix of A consisting of
	 * the first i rows of A */

	fprintf(stderr, "U_upper (%i x %i):\n", U_upper.rows, U_upper.cols);
	matrix_dump(&U_upper, stderr);

	rc = rq_decoder_rfc6330_phase3(rq, A, X, &i, &u);
	test_assert(rc == 0, "rq_decoder_rfc6330_phase2() returned %i", rc);

	/* TODO After this operation, the submatrix of A consisting of the
	 * intersection of the first i rows and columns equals to X, whereas the
	 * matrix U_upper is transformed to a sparse form. */

	/* FIXME - temp - lets just solve this by brute force, then optimize after */


	matrix_free(&U_upper);

	fprintf(stderr, "A (%i x %i):\n", A->rows, A->cols);
	matrix_dump(A, stderr);

	test_assert(matrix_is_identity(A), "A is identity matrix");

	return rc;
}

static int phase_2(rq_t *rq, matrix_t *A, matrix_t *X, int i, int u)
{
	int rc = 0;

	/* create submatricies */
	matrix_t U_upper = matrix_submatrix(A, 0, i, i, u);
	matrix_t U_lower = matrix_submatrix(A, i, i, A->rows - i, u);
	matrix_t I_u = matrix_submatrix(A, i, i, u, u);

#if TEST_DEBUG
	fprintf(stderr, "U (upper):\n");
	matrix_dump(&U_upper, stderr);

	fprintf(stderr, "U (lower):\n");
	matrix_dump(&U_lower, stderr);

	fprintf(stderr, "I_u:\n");
	matrix_dump(&I_u, stderr);
#endif

	rc = rq_decoder_rfc6330_phase2(rq, A, X, &i, &u);
	test_assert(rc == 0, "rq_decoder_rfc6330_phase2() returned %i", rc);

#if TEST_DEBUG
	fprintf(stderr, "A (%i x %i):\n", A->rows, A->cols);
	matrix_dump(A, stderr);

	fprintf(stderr, "X (i x i):\n");
	matrix_dump(X, stderr);
#endif
	if (rc == 0) {
		/* Post Phase-2 tests */
		test_assert(X->rows == i, "X->rows = i");
		test_assert(X->cols == i, "X->cols = i");
		test_assert(A->rows == rq->L, "A->rows = L");
		test_assert(A->cols == rq->L, "A->cols = L");
		test_assert(matrix_is_identity(&I_u), "I_u is identity matrix");
	}
	matrix_free(&U_upper);
	matrix_free(&U_lower);
	matrix_free(&I_u);

	return rc;
}

static int phase_1(rq_t *rq, matrix_t *A, matrix_t *X, int *i, int *u,
		uint8_t *enc, uint32_t ESI[], uint32_t nesi)
{
	uint8_t *dec;
	int rc = 0;

	dec = calloc(rq->K, rq->T);

	/* prepare Matrix A */
	rq_decoder_rfc6330_phase0(rq, A, dec, enc, ESI, nesi);

#if TEST_DEBUG
	fprintf(stderr, "Matrix A (%i x %i)\n", A->rows, A->cols);
	matrix_dump(A, stderr);
#endif
	*i = 0, *u = rq->P;
	rc = rq_decoder_rfc6330_phase1(rq, X, A, i, u);

	test_assert(rc == 0, "rq_decoder_rfc6330_phase1 returned 0");

	fprintf(stderr, "Schedule after Phase 1:\n");
	matrix_schedule_dump(rq->sched, stderr);

	/* Tests at the end of the First Phase: */
	test_assert(X->rows == A->rows, "X.rows == A.rows");
	test_assert(X->cols == A->cols, "X.cols == A.cols");

#if TEST_DEBUG
	matrix_dump(A, stderr);
	test_log("i = %i, u = %i\n", *i, *u);
#endif

	/* The submatrix I defined by the intersection of the first i rows
	 * and first i columns.  This is the identity matrix at the end of each
	 * step in the phase. */
	matrix_t I = matrix_submatrix(A, 0, 0, *i, *i);
	test_assert(matrix_is_identity(&I), "matrix I is identity matrix");
#if TEST_DEBUG
	matrix_dump(&I, stderr);
	assert(I.size == 0);
	matrix_free(&I);
#endif

	/* The submatrix defined by the intersection of the first i rows
	 * and all but the first i columns and last u columns.  All entries of
	 * this submatrix are zero. */
	matrix_t Z0 = matrix_submatrix(A, 0, *i, *i, A->cols - *u - *i);
	test_assert(matrix_is_zero(&Z0), "matrix Z0 is zero");
	matrix_free(&Z0);

	/* The submatrix defined by the intersection of the first i columns
	 * and all but the first i rows.  All entries of this submatrix are
	 * zero. */
	matrix_t Z1 = matrix_submatrix(A, *i, 0, A->rows - *i, *i);
	test_assert(matrix_is_zero(&Z1), "matrix Z1 is zero");
	matrix_free(&Z1);

	/* The phase ends successfully when i + u = L, i.e., when V and the
	 * all zeros submatrix above V have disappeared, and A consists of I,
	 * the all zeros submatrix below I, and U.  The phase ends
	 * unsuccessfully in decoding failure if at some step before V
	 * disappears there is no nonzero row in V to choose in that step. */
	test_assert(*i + *u == rq->L, "i + u == L");

#if TEST_DEBUG
	fprintf(stderr, "X:\n");
	matrix_dump(X, stderr);
#endif
	/* X is supposed to be "lower triangular throughout the first phase",
	 * but it is a straight copy of A which is NOT lower triangular. Only
	 * the i x i portion of X is lower triangular. */
	matrix_t Xii = matrix_submatrix(X, 0, 0, *i, *i);
	test_assert(matrix_is_lower(&Xii), "X (trimmed to i x i) is lower triangular");

	free(dec);

	return rc;
}

static uint8_t *encoder_generate_symbols(rq_t *rq, uint32_t ESI[], int nesi)
{
	uint8_t *enc;
	uint8_t sbn = 0; /* FIXME - hardcoded SBN */
	rq_sym_t sym = {0};

	assert(rq->Z == 1); /* FIXME - only one block supported by this test */
	enc = malloc(nesi * rq->T);

	/* generate random repair symbols */
	sym.sym = enc;
	for (int i = 0; i < nesi; i++) {
		rq_symbol_random(rq, &sym, sbn); /* repair symbols, random order */
		ESI[i] = sym.ESI;
		sym.sym += rq->T;
	}
	return enc;
}

static int encoder_sizetest(uint8_t *src, size_t F, uint16_t T)
{
	rq_t *rq;
	matrix_t C;
	int nesi;
	int rc;

	rq = rq_init(F, T);
	rq_dump(rq, stderr);

	if (F == 0 || T == 0) {
		test_assert(!rq, "zero F or T, rq_init() returns  NULL");
		return 0;
	}
	rc = rq_encode_data(rq, src, F);
	test_assert(rc == 0, "rq_encode_data");

	C = rq_matrix_C_by_SBN(rq, 0);

	nesi = rq->KP + overhead; /* NB: overhead is over and above K', not K */
	if (!rc) {
		uint32_t ESI[nesi];
		uint8_t *enc = encoder_generate_symbols(rq, ESI, nesi);
		matrix_t A, X;
		int i, u;
		rq_t *rq = rq_init(F, T); assert(rq);
		rq->sched = malloc(sizeof(matrix_sched_t));
		memset(rq->sched, 0, sizeof(matrix_sched_t));
		if (!rq->sched) {
			rq_free(rq);
			return -1;
		}

		assert(enc);

		rc = phase_1(rq, &A, &X, &i, &u, enc, ESI, nesi);
		test_assert(rc == 0, "Phase 1");
		rc = phase_2(rq, &A, &X, i, u);
		test_assert(rc == 0, "Phase 2");
		rc = phase_3(rq, &A, &X, i, u);
		test_assert(rc == 0, "Phase 3");
		rc = decodeC(rq, enc, C.base);
		test_assert(rc == 0, "Decode C");
		matrix_free(&X);
		matrix_free(&A);
		rq_free(rq);
		free(enc);
	}
	rq_free(rq);
	return rc;
}

int main(void)
{
	uint8_t *srcobj = NULL;
	int rc;

	loginit();
	test_name("5.4.2.2 RFC Decoder");
	for (overhead = OMIN; overhead <= OMAX; overhead++) {
		for (int T = TMIN; T <= TMAX; T *= RQ_AL) {
			for (size_t F = FMIN; F <= FMAX; F++) {
				int ok = 0; int tests = 0;
				for (int i = 0; i < REPS; i++) {
					srcobj = generate_source_object(F);
					rc = encoder_sizetest(srcobj, F, T);
					free(srcobj);
					if (!rc) {
						test_log("SUCCESS with F=%zu, T=%u\n", F, T);
						ok++;
					}
					tests++;
					//if (rc) goto test_done; /* stop on first failure */
				}
				test_assert(ok == tests, "F=%zu, T=%u %i/%i (overhead = %u) tests ok",
						F, T, ok, tests, overhead);
			}
		}
	}
//test_done:
	test_log("test done\n");

	return fails;
}
