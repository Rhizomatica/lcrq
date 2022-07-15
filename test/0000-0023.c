/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <lcrq.h>

int main(void)
{
	rq_pid_t pid = 0;
	uint32_t esi = 19;
	uint8_t sbn = 42;

	loginit();
	test_name("PID/SBN/ESI bitshifting macros");

	for (int i = 0; i < 255; i++)
	{
		sbn = test_randomnumber(1 << 8);
		esi = test_randomnumber(1 << 24);
		pid = rq_pidsetsbn(pid, sbn);
		pid = rq_pidsetesi(pid, esi);
		test_assert(sbn == rq_pid2sbn(pid), "sbn");
		test_assert(esi == rq_pid2esi(pid), "esi");
	}

	return fails;
}
