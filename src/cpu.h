/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022-2023 Brett Sheffield <bacs@librecast.net> */
/*
 * cpu.h - test CPU instruction set support
 *
 * Acknowledgements:
 *
 *  Agner, the BLAKE3 team and Intel, whose code and manuals
 *  were all used as references.
 */
#ifndef _CPU_H
#define _CPU_H 1

#include "config.h"
#ifdef HAVE_IMMINTRIN_H
# include <immintrin.h>
#endif
#include <stdint.h>

#define LCRQ_NOSIMD		0
#define LCRQ_SIMD_MMX		1
#define LCRQ_SIMD_SSE		2
#define LCRQ_SIMD_SSE2		3
#define LCRQ_SIMD_SSE3		4
#define LCRQ_SIMD_SSSE3		5
#define LCRQ_SIMD_SSE4_1	6
#define LCRQ_SIMD_SSE4_2	7
#define LCRQ_SIMD_POPCNT	8
#define LCRQ_SIMD_AVX		9
#define LCRQ_SIMD_AVX2		10

#define MMX    1 << LCRQ_SIMD_MMX
#define SSE    1 << LCRQ_SIMD_SSE
#define SSE2   1 << LCRQ_SIMD_SSE2
#define SSE3   1 << LCRQ_SIMD_SSE3
#define SSSE3  1 << LCRQ_SIMD_SSSE3
#define SSE4_1 1 << LCRQ_SIMD_SSE4_1
#define SSE4_2 1 << LCRQ_SIMD_SSE4_2
#define POPCNT 1 << LCRQ_SIMD_POPCNT
#define AVX    1 << LCRQ_SIMD_AVX
#define AVX2   1 << LCRQ_SIMD_AVX2

#ifdef HAVE_IMMINTRIN_H
#if defined ( NOSIMD_COMMON )
#elif defined ( __AVX2__ )
# define LCRQ_SIMD LCRQ_SIMD_AVX2
#elif defined ( __AVX__ )
# define LCRQ_SIMD LCRQ_SIMD_AVX
#elif defined ( __POPCNT__ )
# define LCRQ_SIMD LCRQ_SIMD_POPCNT
#elif defined ( __SSE4_2__ )
# define LCRQ_SIMD LCRQ_SIMD_SSE4_2
#elif defined ( __SSE4_1__ )
# define LCRQ_SIMD LCRQ_SIMD_SSE4_1
#elif defined ( __SSSE3__ )
# define LCRQ_SIMD LCRQ_SIMD_SSSE3
#elif defined ( __SSE3__ )
# define LCRQ_SIMD LCRQ_SIMD_SSE3
#elif defined ( __SSE2__ )
# define LCRQ_SIMD LCRQ_SIMD_SSE2
#elif defined ( __SSE__ )
# define LCRQ_SIMD LCRQ_SIMD_SSE
#elif defined ( __MMX__ )
# define LCRQ_SIMD LCRQ_SIMD_MMX
#endif
#endif

int cpu_instruction_set(void);
#endif
