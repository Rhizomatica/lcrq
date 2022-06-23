/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include <matrix.h>
#include <gf256.h>
#include <stdlib.h>
#include <string.h>
#include <sys/param.h>
#include <unistd.h>

#ifdef INTEL_SSE3
#include <emmintrin.h>
#include <immintrin.h>
#endif

#define VSZ 256

uint8_t reclen[5] = {
	0, /* NOOP */
	sizeof(matrix_op_swap_t),
	sizeof(matrix_op_swap_t),
	sizeof(matrix_op_add_t),
	sizeof(matrix_op_mul_t)
};

matrix_t *matrix_new(matrix_t *m, const int rows, const int cols, uint8_t *base, int flags)
{
	m->rows = rows;
	m->cols = cols;
	m->cvec = ((flags & MATRIX_VEC) == MATRIX_VEC) ? roundup(cols, VSZ) : 0;
	m->trans = 0;
	m->roff = 0;
	m->coff = 0;
	if (base || !m->cvec) {
		m->stride = (size_t)cols * sizeof(uint8_t);
		m->size = (size_t)rows * (size_t)cols * sizeof(uint8_t);
	}
	else {
		m->stride = m->cvec * sizeof(uint8_t);
		m->size = (size_t)rows * m->cvec * sizeof(uint8_t);
	}
	m->base = (base) ? base : malloc(m->size);
	return m;
}

matrix_t matrix_submatrix(const matrix_t *A, const int off_rows, const int off_cols,
		const int rows, const int cols)
{
	const long int r = (A->trans) ? cols : rows;
	const long int c = (A->trans) ? rows : cols;
	const long int roff = (A->trans) ? off_cols : off_rows;
	const long int coff = (A->trans) ? off_rows : off_cols;
	assert(rows <= matrix_rows(A));
	assert(cols <= matrix_cols(A));
	matrix_t sub = {
		.rows = r,
		.cols = c,
		.base = A->base + roff * A->stride + coff,
		.trans = A->trans,
		.stride = A->stride,
		.roff = roff,
		.coff = coff,
	};
	return sub;
}

uint8_t *matrix_zero_row(matrix_t *m, int row)
{
	return memset(matrix_ptr_row(m, row), 0, m->cvec);
}

/* A submatrix needs to be zeroed row by row.
 * A full matrix has a size which can be passed to memset */
void matrix_zero(matrix_t *m)
{
	assert(m && m->base);
	if (m->size) memset(m->base, 0, m->size);
	else for (int i = 0; i < m->rows; i++) matrix_zero_row(m, i);
}

matrix_t *matrix_identity(matrix_t *m)
{
	assert(m->rows == m->cols);
	matrix_zero(m);
	for (int i = 0; i < m->rows; i++) {
		matrix_set(m, i, i, 1);
	}
	return m;
}

/* return nonzero if m is an identity matrix */
int matrix_is_identity(matrix_t *m)
{
	for (int i = 0; i < m->rows; i++) {
		for (int j = 0; j < m->cols; j++) {
			const uint8_t x = matrix_get(m, i, j);
			const uint8_t y = (i == j) ? 1 : 0;
			if (x != y) return 0;
		}
	}
	return -1;
}

/* return nonzero if all the entries above the main diagonal are zero */
int matrix_is_lower(matrix_t *m)
{
	for (int i = 0; i < m->rows; i++) {
		for (int j = i + 1; j < m->cols; j++) {
			if (matrix_get(m, i, j)) return 0;
		}
	}
	return -1;
}

int matrix_is_zero(matrix_t *m)
{
	for (int i = 0; i < m->rows; i++) {
		for (int j = 0; j < m->cols; j++) {
			if (matrix_get_s(m, i, j)) return 0;
		}
	}
	return -1;
}

int matrix_row_degree(matrix_t *m, int row)
{
	int d = 0;
	for (int j = 0; j < matrix_cols(m); j++) {
		if (matrix_get(m, row, j)) d++;
	}
	return d;
}

void matrix_dump(matrix_t *mat, FILE *stream)
{
	fprintf(stream, "\n");
	for (int r = 0; r < matrix_rows(mat); r++) {
		for (int c = 0; c < matrix_cols(mat); c++) {
			const uint8_t v = matrix_get(mat, r, c);
			if (!v) fprintf(stream, " --");
			else fprintf(stream, " %02x", matrix_get(mat, r, c));
		}
		fprintf(stream, "\n");
	}
	fprintf(stream, "\n");
}

void matrix_inc_gf256(matrix_t *mat, const int row, const int col, const uint8_t val)
{
	const uint8_t sum = matrix_get(mat, row, col) ^ val;
	matrix_set(mat, row, col, sum);
}

matrix_t *matrix_multiply_inplace(const matrix_t *x, matrix_t *y)
{
	const int yrows = matrix_rows(y);
	const int xcols = matrix_cols(x);

	assert(xcols == yrows);

	for (int i = 0; i < yrows; i++) {
		for (int j = 0; j < xcols; j++) {
			for (int k = 0; k < xcols; k++) {
				uint8_t *a = MADDR(y, k, j);
				const uint8_t b = matrix_get(x, i, k);
				*a ^= GF256MUL(*a, b);
			}
		}
	}
	return y;
}

/* GF(256) dot product of x and y returned in p. Allocate p->base if required */
matrix_t *matrix_multiply_gf256(const matrix_t *x, const matrix_t *y, matrix_t *p)
{
	const int xrows = matrix_rows(x);
	const int yrows = matrix_rows(y);
	const int xcols = matrix_cols(x);
	const int ycols = matrix_cols(y);
	assert(xcols == yrows);
	if (!p->base) {
		matrix_new(p, xrows, ycols, NULL, 0);
	}
	else {
		if (p->cols != xrows || p->rows != yrows)
			return NULL;
	}
	matrix_zero(p);
	for (int i = 0; i < p->rows; i++) {
		for (int j = 0; j < p->cols; j++) {
			uint8_t v = 0;
			for (int k = 0; k < xcols; k++) {
				const uint8_t a = matrix_get(x, i, k);
				const uint8_t b = matrix_get(y, k, j);
				v ^= GF256MUL(a, b);
			}
			matrix_set(p, i, j, v);
		}
	}

	return p;
}

matrix_t *matrix_swap_rows(matrix_t *m, const int r1, const int r2)
{
	for (int i = 0; i < m->cols; i++) {
		uint8_t *a = MADDR(m, r1, i);
		uint8_t *b = MADDR(m, r2, i);
		SWAP(*a, *b);
	}
	return m;
}

matrix_t *matrix_swap_cols(matrix_t *m, const int c1, const int c2)
{
	for (int i = 0; i < m->rows; i++) {
		uint8_t *a = MADDR(m, i, c1);
		uint8_t *b = MADDR(m, i, c2);
		SWAP(*a, *b);
	}
	return m;
}

void matrix_row_add(matrix_t *dst, const int drow, const matrix_t *src, const int srow)
{
	assert(matrix_cols(dst) == matrix_cols(src));
	uint8_t *d = matrix_ptr_row(dst, drow);
	uint8_t *s = matrix_ptr_row(src, srow);
	const int mcols = matrix_cols(dst);
#ifdef INTEL_SSE3
	const int mod = mcols % 16;
	const int maxv = mcols - mod;
	int j;
	for (j = 0; j < maxv; j += 16) {
		__m128i S = _mm_loadu_si128((const __m128i_u *)&s[j]);
		__m128i D = _mm_loadu_si128((const __m128i_u *)&d[j]);
		D = _mm_xor_si128(D, S);
		_mm_storeu_si128((__m128i*)&d[j], D);
	}
#else
	int j = 0;
#endif
	for (; j < mcols; j++) d[j] ^= s[j];
}

matrix_t matrix_add(const matrix_t *x, const matrix_t *y)
{
	assert(matrix_rows(x) == matrix_rows(y));
	assert(matrix_cols(x) == matrix_cols(y));
	matrix_t p = matrix_dup(x);
	matrix_zero(&p);
	for (int i = 0; i < matrix_rows(x); i++) {
		matrix_row_add(&p, i, x, i);
		matrix_row_add(&p, i, y, i);
	}
	return p;
}

void matrix_row_add_val(matrix_t *m, const int row, const uint8_t val)
{
	for (int col = 0; col < m->cols; col++) {
		matrix_set(m, row, col, matrix_get(m, row, col) ^ val);
	}
}

void matrix_row_div(matrix_t *m, const int row, const uint8_t val)
{
	for (int col = 0; col < m->cols; col++) {
		uint8_t a = matrix_get(m, row, col);
		uint8_t b = val;
		uint8_t v = GF256DIV(a, b);
		matrix_set(m, row, col, v);
	}
}

void matrix_col_mul(matrix_t *m, const int col, const int off, const uint8_t v)
{
	for (int row = off; row < matrix_rows(m); row++) {
		matrix_set(m, row, col, GF256MUL(matrix_get(m, row, col), v));
	}
}

void matrix_col_copy(matrix_t *dst, const int dcol, const matrix_t *src, const int scol)
{
	for (int row = 0; row < src->rows; row++) {
		matrix_set(dst, row, dcol, matrix_get(src, row, scol));
	}
}

uint8_t * matrix_row_copy(matrix_t *dst, const int drow, const matrix_t *src, const int srow)
{
	uint8_t *dptr, *sptr;
	dptr = matrix_ptr_row(dst, drow);
	sptr = matrix_ptr_row(src, srow);
	return memcpy(dptr, sptr, src->stride);
}

int matrix_pivot(matrix_t *A, int j, int P[], int Q[])
{
	const int Arows = matrix_rows(A);
	const int Acols = matrix_cols(A);
	for (int col = j; col < Acols; col++) {
		for (int row = j; row < Arows; row++) {
			if (matrix_get_s(A, row, j)) {
				/* pivot found, move in place, update P+Q */
				if (row != j) {
					matrix_swap_rows(A, row, j);
					SWAP(P[row], P[j]);
				}
				if (col != j) {
					matrix_swap_cols(A, col, j);
					SWAP(Q[col], Q[j]);
				}
				return 1;
			}
		}
	}
	return 0;
}

int matrix_LU_decompose(matrix_t *A, int P[], int Q[])
{
	const int Arows = matrix_rows(A);
	const int Acols = matrix_cols(A);
	const int n = MIN(Arows, Acols);
	int i;

	/* initialize permutations matricies */
	for (i = 0; i < Arows; i++) P[i] = i;
	for (i = 0; i < Acols; i++) Q[i] = i;

	/* LU decomposition */
	for (i = 0; i < n; i++) {
		uint8_t *Aii = MADDR(A, i, i);
		if (!matrix_pivot(A, i, P, Q)) break;
		for (int j = i + 1; j < Arows; j++) {
			uint8_t *Aji = MADDR(A, j, i);
			*Aji = GF256DIV(*Aji, *Aii);
			if (*Aji) for (int k = i + 1; k < Acols; k++) {
				const uint8_t b = matrix_get_s(A, i, k);
				A->base[k + j * A->stride] ^= GF256MUL(*Aji, b);
			}
		}
	}
	return i;
}

void matrix_inverse_LU(matrix_t *IA, const matrix_t *LU, const int P[])
{
	int n = MIN(matrix_rows(LU), matrix_cols(LU));

	assert(n == matrix_cols(LU));

	if (!IA->base) matrix_new(IA, matrix_rows(LU), matrix_cols(LU), NULL, 0);

	for (int j = 0; j < matrix_cols(IA); j++) {
		for (int i = 0; i < matrix_cols(IA); i++) {
			uint8_t v = (P[i] == j) ? 1 : 0;
			matrix_set(IA, i, j, v);
			for (int k = 0; k < i; k++) {
				const uint8_t a = matrix_get(LU, i, k);
				const uint8_t b = matrix_get(IA, k, j);
				const uint8_t ab = GF256MUL(a, b);
				matrix_inc_gf256(IA, i, j, ab);
			}
		}

		for (int i = n - 1; i >= 0; i--) {
			for (int k = i + 1; k < matrix_cols(IA); k++) {
				const uint8_t a = matrix_get(LU, i, k);
				const uint8_t b = matrix_get(IA, k, j);
				const uint8_t ab = GF256MUL(a, b);
				matrix_inc_gf256(IA, i, j, ab);
			}
			const uint8_t v = gf256_div(matrix_get(IA, i, j), matrix_get(LU, i, i));
			matrix_set(IA, i, j, v);
		}
	}
}

void matrix_row_mul(matrix_t *m, const int row, const int off, const uint8_t val)
{
	uint8_t *d = matrix_ptr_row(m, row) + off;
	for (int col = 0; col < m->cols; col++) {
		d[col] = GF256MUL(d[col], val);
	}
}

#ifdef INTEL_SSE3
inline static __m128i mul_128(__m128i A, uint8_t y)
{
	__m128i table1 = _mm_loadu_si128((const __m128i_u *)GF256LR[y][0]);
	__m128i table2 = _mm_loadu_si128((const __m128i_u *)GF256LR[y][1]);
	__m128i mask1 = _mm_set1_epi8((uint8_t)0x0f);
	__m128i mask2 = _mm_set1_epi8((uint8_t)0xf0);
	__m128i l, h;

	l = _mm_and_si128(A, mask1);
	l = _mm_shuffle_epi8(table1, l);
	h = _mm_and_si128(A, mask2);
	h = _mm_srli_epi64(h, 4);
	h = _mm_shuffle_epi8(table2, h);
	return _mm_xor_si128(h, l);
}
#endif

void matrix_row_mul_byrow(matrix_t *m, const int rdst, const int off, const int rsrc, const uint8_t y)
{
	assert(y);
	uint8_t *d = matrix_ptr_row(m, rdst) + off;
	uint8_t *s = matrix_ptr_row(m, rsrc) + off;
#ifdef INTEL_SSE3
	const int max = m->cols - off;
	const int mod = max % 16;
	const int maxv = max - mod;
	int i;

	for (i = 0; i < maxv; i += 16) {
		__m128i S = _mm_loadu_si128((const __m128i_u *)&s[i]);
		__m128i D = _mm_loadu_si128((const __m128i_u *)&d[i]);
		S = mul_128(S, y);
		D = _mm_xor_si128(D, S);
		_mm_storeu_si128((__m128i*)&d[i], D);
	}
#else
	const int max = m->cols - off;
	int i = 0;
#endif
	for (; i < max; i++) {
		d[i] ^= GF256MUL(s[i], y);
	}
}

void matrix_solve_LU(matrix_t *X, const matrix_t *Y, const matrix_t *LU, const int P[], const int Q[])
{
	int n = MIN(matrix_rows(LU), matrix_cols(LU));

	if (!X->base) matrix_new(X, matrix_rows(LU), matrix_cols(Y), NULL, 0);

	assert(matrix_cols(LU) == matrix_rows(X));
	assert(matrix_rows(LU) == matrix_rows(Y));
	assert(matrix_cols(X) == matrix_cols(Y));

	if (!X->base) matrix_new(X, matrix_rows(LU), matrix_cols(LU), NULL, 0);

	for (int i = 0; i < n; i++) {
		matrix_row_copy(X, Q[i], Y, P[i]);
		for (int j = 0; j < i; j++) {
			const uint8_t LUij = matrix_get_s(LU, i, j);
			if (LUij) matrix_row_mul_byrow(X, Q[i], 0, Q[j], LUij);
		}
	}
	for (int i = n - 1; i >= 0; i--) {
		for (int j = i + 1; j < matrix_cols(LU); j++) {
			const uint8_t LUij = matrix_get_s(LU, i, j);
			if (LUij) matrix_row_mul_byrow(X, Q[i], 0, Q[j], LUij);
		}
		matrix_row_mul(X, Q[i], 0, GF256INV(matrix_get_s(LU, i, i)));
	}
}

static int matrix_pivot_sched(matrix_t *A, int j, matrix_sched_t *sched)
{
	const int Arows = matrix_rows(A);
	const int Acols = matrix_cols(A);
	for (int col = j; col < Acols; col++) {
		for (int row = j; row < Arows; row++) {
			if (matrix_get_s(A, row, j)) {
				/* pivot found, move in place, update P+Q */
				if (row != j) {
					matrix_swap_rows(A, row, j);
					matrix_sched_row(sched, A->roff + row, A->roff + j);
				}
				if (col != j) {
					matrix_swap_cols(A, col, j);
					matrix_sched_row(sched, A->coff + col, A->coff + j);
				}
				return 1;
			}
		}
	}
	return 0;
}

/* reduce matrix to identity, track operations in matrix schedule, return rank */
int matrix_gauss_elim(matrix_t *A, matrix_sched_t *sched)
{
	int j;

	/* lower echelon form */
	for (j = 0; j < A->cols; j++) {
		if (!matrix_pivot_sched(A, j, sched)) break;
		/* first, reduce the pivot row so jj = 1 */
		uint8_t jj = matrix_get_s(A, j, j);
		if (jj != 1) {
			const uint8_t b = GF256INV(jj);
			if (b) {
				matrix_row_mul(A, j, 0, b);
				matrix_sched_mul(sched, A->roff + j, b);
			}
		}
		for (int i = j + 1; i < A->rows; i++) {
			/* add pivot row (j) * factor to row i so that ij == 0 */
			jj = matrix_get_s(A, j, j);
			const uint8_t ij = matrix_get_s(A, i, j);
			if (ij) {
				const uint8_t f = gf256_div(ij, jj);
				if (f) {
					matrix_row_mul_byrow(A, i, 0, j, f);
					matrix_sched_add(sched, A->roff + i, 0, A->roff + j, f);
				}
			}
		}
	}
	/* finish upper triangle */
	for (int i = 0; i < A->cols; i++) {
		for (int j = i + 1; j < A->cols; j++) {
			const uint8_t ij = matrix_get_s(A, i, j);
			if (ij) {
				const uint8_t jj = matrix_get_s(A, j, j);
				const uint8_t f = gf256_div(ij, jj);
				if (f) {
					matrix_row_mul_byrow(A, i, 0, j, f);
					matrix_sched_add(sched, A->roff + i, 0, A->roff + j, f);
				}
			}
		}
	}

	return j; /* rank */
}

matrix_t *matrix_inverse(matrix_t *A, matrix_t *I)
{
	matrix_new(I, A->rows, A->cols, NULL, 0);
	matrix_identity(I);

	/* lower echelon form */
	for (int j = 0; j < A->cols; j++) {
		/* first, reduce the pivot row so jj = 1 */
		uint8_t jj = matrix_get(A, j, j);
		assert(jj); /* without row swaps, cannot proceed */
		if (jj != 1) {
			matrix_row_div(A, j, jj);
			matrix_row_div(I, j, jj);
		}
		for (int i = j + 1; i < A->rows; i++) {
			/* add pivot row (j) * factor to row i so that ij == 0 */
			jj = matrix_get(A, j, j);
			const uint8_t ij = matrix_get(A, i, j);
			const uint8_t f = gf256_div(ij, jj);
			matrix_row_mul_byrow(A, i, 0, j, f);
			matrix_row_mul_byrow(I, i, 0, j, f);
		}
	}
	/* finish upper triangle */
	for (int i = 0; i < A->cols; i++) {
		for (int j = i + 1; j < A->cols; j++) {
			const uint8_t ij = matrix_get(A, i, j);
			const uint8_t jj = matrix_get(A, j, j);
			const uint8_t f = gf256_div(ij, jj);
			matrix_row_mul_byrow(A, i, 0, j, f);
			matrix_row_mul_byrow(I, i, 0, j, f);
		}
	}

	return I;
}

matrix_t *matrix_copy(matrix_t *dst, const matrix_t *src)
{
	if (src->size) {
		memcpy(dst->base, src->base, src->size);
		dst->size = src->size;
	}
	else {
		/* size = 0 => submatrix, copy row by row */
		for (int i = 0; i < matrix_rows(src); i++) {
			matrix_row_copy(dst, i, src, i);
		}
	}
	return dst;
}

matrix_t matrix_dup(const matrix_t *src)
{
	matrix_t m = {0};
	matrix_new(&m, matrix_rows(src), matrix_cols(src), NULL, 0);
	matrix_copy(&m, src);
	return m;
}

void matrix_transpose(matrix_t *mat)
{
	mat->trans = !(mat->trans);
}

static matrix_op_t *matrix_sched_op(matrix_sched_t *sched, uint8_t optype)
{
	uint8_t lasttype;
	matrix_op_t *op;
resized:
	lasttype = *(sched->last) & 0x0f;
	op = (matrix_op_t *)(sched->last + reclen[lasttype]);
	if (sched->len <= (size_t)(sched->last - sched->base + reclen[lasttype] + reclen[optype])) {
		matrix_schedule_resize(sched);
		goto resized;
	}
	op->swp.type = (lasttype << 4) | optype;
	sched->last = (uint8_t *)op;
	(sched->ops)++;
	return op;
}

void matrix_sched_add(matrix_sched_t *sched, uint16_t dst, uint16_t off, uint16_t src, uint8_t beta)
{
	matrix_op_add_t *op = (matrix_op_add_t *)matrix_sched_op(sched, MATRIX_OP_ADD);
	op->dst = dst;
	op->src = src;
	op->beta = beta;
	op->off = off;
}

void matrix_sched_mul(matrix_sched_t *sched, uint16_t dst, uint8_t beta)
{
	matrix_op_mul_t *op = (matrix_op_mul_t *)matrix_sched_op(sched, MATRIX_OP_MUL);
	op->dst = dst;
	op->beta = beta;
}

static void matrix_sched_swap(matrix_sched_t *sched, uint16_t a, uint16_t b, uint8_t optype)
{
	matrix_op_swap_t *op = (matrix_op_swap_t *)matrix_sched_op(sched, optype);
	op->a = a;
	op->b = b;
}

void matrix_sched_row(matrix_sched_t *sched, uint16_t a, uint16_t b)
{
	matrix_sched_swap(sched, a, b, MATRIX_OP_ROW);
}

void matrix_sched_col(matrix_sched_t *sched, uint16_t a, uint16_t b)
{
	matrix_sched_swap(sched, a, b, MATRIX_OP_COL);
}

void matrix_schedule_dump(matrix_sched_t *sched, FILE *stream)
{
	uint8_t type;
	matrix_op_t *o;
	fprintf(stream, "-- schedule begins --\n");
	for (uint8_t *op = sched->base; *op; op += reclen[type]) {
		o = (matrix_op_t *)op;
		type = (*op) & 0x0f;
		switch (type) {
		case MATRIX_OP_NOOP:
			fprintf(stream, "NOOP\n");
			break;
		case MATRIX_OP_ROW:
			fprintf(stream, "OP_ROW %u %u\n", o->swp.a, o->swp.b);
			break;
		case MATRIX_OP_COL:
			fprintf(stream, "OP_COL %u %u\n", o->swp.a, o->swp.b);
			break;
		case MATRIX_OP_ADD:
			fprintf(stream, "OP_ADD %u %u %u %u\n", o->add.dst, o->add.off,
					o->add.src, o->add.beta);
			break;
		case MATRIX_OP_MUL:
			fprintf(stream, "OP_MUL %u %u\n", o->mul.dst, o->mul.beta);
			break;
		default:
			assert(0); /* MUST not occur */
		}
	}
	fprintf(stream, "-- schedule ends -- (%zu bytes, %zu ops)\n", sched->len, sched->ops);
}

void matrix_schedule_reverse(matrix_t *m, matrix_sched_t *sched)
{
	matrix_op_t *o;
	uint8_t *op = sched->last;
	uint8_t type;

	fprintf(stderr, "-- reverse replay begins --\n");
	type = (*op) & 0x0f;
	do {
		o = (matrix_op_t *)op;
		switch (type) {
		case MATRIX_OP_NOOP:
			fprintf(stderr, "NOOP\n");
			break;
		case MATRIX_OP_ROW:
			fprintf(stderr, "OP_ROW %u %u\n", o->swp.a, o->swp.b);
			matrix_swap_rows(m, o->swp.a, o->swp.b);
			break;
		case MATRIX_OP_COL:
			fprintf(stderr, "OP_COL %u %u\n", o->swp.a, o->swp.b);
			matrix_swap_cols(m, o->swp.a, o->swp.b);
			break;
		case MATRIX_OP_ADD:
			fprintf(stderr, "OP_ADD %u %u %u %u\n", o->add.dst, o->add.off,
					o->add.src, o->add.beta);
			matrix_row_mul_byrow(m, o->add.dst, o->add.off, o->add.src, o->add.beta);
			break;
		case MATRIX_OP_MUL:
			fprintf(stderr, "OP_MUL %u %u\n", o->mul.dst, o->mul.beta);
			matrix_row_mul(m, o->mul.dst, 0, o->mul.beta);
			break;
		default:
			assert(0); /* MUST not occur */
		}
		type = (*op) >> 4; /* type of previous record */
		op -= reclen[type];
	} while (type);
	fprintf(stderr, "-- replay ends -- (%zu bytes, %zu ops)\n", sched->len, sched->ops);
}

void matrix_schedule_replay(matrix_t *m, matrix_sched_t *sched)
{
	uint8_t type;
	matrix_op_t *o;
	fprintf(stderr, "-- replay begins --\n");
	for (uint8_t *op = sched->base; *op; op += reclen[type]) {
		o = (matrix_op_t *)op;
		type = (*op) & 0x0f;
		switch (type) {
		case MATRIX_OP_NOOP:
			fprintf(stderr, "NOOP\n");
			break;
		case MATRIX_OP_ROW:
			fprintf(stderr, "OP_ROW %u %u\n", o->swp.a, o->swp.b);
			matrix_swap_rows(m, o->swp.a, o->swp.b);
			break;
		case MATRIX_OP_COL:
			fprintf(stderr, "OP_COL %u %u\n", o->swp.a, o->swp.b);
			matrix_swap_cols(m, o->swp.a, o->swp.b);
			break;
		case MATRIX_OP_ADD:
			fprintf(stderr, "OP_ADD %u %u %u %u\n", o->add.dst, o->add.off,
					o->add.src, o->add.beta);
			matrix_row_mul_byrow(m, o->add.dst, o->add.off, o->add.src, o->add.beta);
			break;
		case MATRIX_OP_MUL:
			fprintf(stderr, "OP_MUL %u %u\n", o->mul.dst, o->mul.beta);
			matrix_row_mul(m, o->mul.dst, 0, o->mul.beta);
			break;
		default:
			assert(0); /* MUST not occur */
		}
	}
	fprintf(stderr, "-- replay ends -- (%zu bytes, %zu ops)\n", sched->len, sched->ops);
}

void matrix_schedule_free(matrix_sched_t *sched)
{
	if (sched) {
		free(sched->base);
		sched->base = NULL;
	}
}

uint8_t *matrix_schedule_resize(matrix_sched_t *sched)
{
	uint8_t *ptr;
	long sz;
	if ((sz = sysconf(_SC_PAGESIZE)) == -1) return NULL;
	ptr = realloc(sched->base, sched->len + (size_t)sz);
	if (ptr) {
		if (sched->base) sched->last = ptr + (sched->last - sched->base);
		else sched->last = ptr;
		sched->base = ptr;
		memset(ptr + sched->len, 0, sz);
		sched->len += sz;
	}
	return ptr;
}

uint8_t *matrix_schedule_init(matrix_sched_t *sched)
{
	sched->base = matrix_schedule_resize(sched);
	if (sched->base) sched->last = sched->base;
	return sched->base;
}

void matrix_free(matrix_t *m)
{
	if (m->size) free(m->base);
	m->base = NULL;
}
