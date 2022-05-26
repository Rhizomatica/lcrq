/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#ifndef GF256_H
#define GF256_H

#include <stdint.h>

/* GF(256) operations as per RFC6330 (5.7) */

/* add (or subtract) in GF(256) */
uint8_t gf256_add(uint8_t a, uint8_t b);

uint8_t gf256_exp(uint8_t e);

uint8_t gf256_log(uint8_t v);

uint8_t gf256_inv(uint8_t v);

uint8_t gf256_mul(uint8_t a, uint8_t b);

#endif /* GF256_H */
