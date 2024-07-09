/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022-2023 Brett Sheffield <bacs@librecast.net> */

#include "cpu.h"

static inline uint64_t xgetbv (int ctr)
{
	uint32_t a, d;
	__asm("xgetbv" : "=a"(a),"=d"(d) : "c"(ctr) : );
	return a | (((uint64_t)d) << 32);
}

static inline void cpuid(uint32_t e_x[4])
{
	__asm__ __volatile__("cpuid" : "+a"(e_x[0]), "+b"(e_x[1]), "+c"(e_x[2]), "=d"(e_x[3]) : );
}

#define eax e_x[0]
#define ebx e_x[1]
#define ecx e_x[2]
#define edx e_x[3]
int cpu_instruction_set(void) {
	uint32_t e_x[4] = {0}; /* registers eax, ebx, ecx, edx */
	static int isets = -1;
	if (isets >= 0) return isets;
	isets = 0;
	cpuid(e_x);
	if (!eax) return isets;
	eax = 1, cpuid(e_x);

	/* check edx register */
	if (!(edx & (1 << 0))) return isets;
	/* floating point */
	if (!(edx & (1 << 23))) return isets;
	/* MMX */
	isets |= MMX;
	if (!(edx & (1 << 15))) return isets;
	/* conditional move */
	if (!(edx & (1 << 24))) return isets;
	/* FXSAVE */
	if (!(edx & (1 << 25))) return isets;
	/* SSE */
	isets |= SSE;
	if (!(edx & (1 << 26))) return isets;
	/* SSE2 */
	isets |= SSE2;

	/* check ecx register */
	if (!(ecx & (1 << 0))) return isets;
	/* SSE3 */
	isets |= SSE3;
	if (!(ecx & (1 << 9))) return isets;
	/* SSSE3 */
	isets |= SSSE3;
	if (!(ecx & (1 << 19))) return isets;
	/* SSE4.1 */
	isets |= SSE4_1;
	if (!(ecx & (1 << 23))) return isets;
	/* POPCNT */
	isets |= POPCNT;
	if (!(ecx & (1 << 20))) return isets;
	/* SSE4.2 */
	isets |= SSE4_2;
	if (!(ecx & (1 << 27))) return isets;
	/* OSXSAVE */
	if ((xgetbv(0) & 6) != 6) return isets;
	if (!(ecx & (1 << 28))) return isets;
	/* AVX */
	isets |= AVX;

	/* check ebx register */
	eax = 7, ecx = 0, cpuid(e_x);
	if (!(ebx & (1 << 5))) return isets;
	/* AVX2 */
	isets |= AVX2;

	return isets;
}
