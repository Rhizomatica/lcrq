/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#ifndef LCRQ_H
#define LCRQ_H 1

#include <assert.h>
#include <stddef.h>
#include <stdint.h>
#ifndef NDEBUG
#include <stdio.h>
#endif

#define RQ_AL 4

/* allow for K + RQ_OVERHEAD repair symbols */
#define RQ_OVERHEAD 5

/* NB: ESI is a 24-bit unsigned integer (3.2) */
#define RQ_ESI_MAX 0xffffff

/* default working memory value */
extern size_t RQ_WS_DEFAULT;

/* rounding integer division */
#define CEIL(x, y)  (((x) + ((y) - 1)) / (y))
#define FLOOR(x, y) ((x) / (y))

enum {
	RQ_SOURCE = 1,
	RQ_REPAIR = 2,
	RQ_REPEAT = 4,
	RQ_ASYNC  = 8,
	RQ_RAND   = 16,
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

/* FEC Payload IDs (3.2):
 * 0                   1                   2                   3
 * 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1
 * +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 * |     SBN       |               Encoding Symbol ID              |
 * +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 */
typedef uint32_t rq_pid_t;

/* FEC Object Transmission Information (3.3):
 * 0                   1                   2                   3
 * 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1
 * +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 * |                      Transfer Length (F)                      |
 * +               +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 * |               |    Reserved   |           Symbol Size (T)     |
 * +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 */
#if 0
typedef uint64_t rq_oti_t;
#endif

/* Scheme-Specific (3.3.3):
 * 0                   1                   2                   3
 * 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1 2 3 4 5 6 7 8 9 0 1
 * +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 * |       Z       |              N                |       Al      |
 * +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
 */
#if 0
typedef uint32_t rq_scheme_t;
#endif

/* rq_init - initialize RaptorQ context
 * creates and returns a new RaptorQ context and sets up the environment.
 * Call rq_free(3) when done.*/
rq_t *rq_init(const uint64_t F, const uint16_t T);

/* rq_free - free RaptorQ context */
void rq_free(rq_t *rq);

/* Query RaptorQ internal values */
uint64_t rq_F(const rq_t * const rq);
uint16_t rq_T(const rq_t * const rq);
uint16_t rq_Z(const rq_t * const rq);
uint16_t rq_N(const rq_t * const rq);
uint8_t rq_Al(const rq_t * const rq);
uint16_t rq_KP(const rq_t * const rq);
uint16_t rq_K(const rq_t * const rq);

int rq_encode(rq_t *rq, uint8_t *data, const size_t len);
int rq_decode(rq_t *rq, uint8_t *dec, uint8_t *enc, uint32_t ESI[], uint32_t nesi);

int __attribute__((deprecated))rq_encode_data_rfc(rq_t *rq, uint8_t *data, const size_t len);
int __attribute__((deprecated))rq_decode_block_rfc(rq_t *rq, uint8_t *dec, uint8_t *enc, uint32_t ESI[], uint32_t nesi);

uint8_t *rq_pkt_gen(const rq_t *rq, rq_pid_t *pid, uint8_t *sym, int flags);

#endif /* LCRQ_PVT_H */
