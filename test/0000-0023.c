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
		sbn = test_randomnumber(UINT8_MAX);
		esi = test_randomnumber(UINT32_MAX >> 8);
		pid = rq_pidset(sbn, esi);
		test_log("%03i: sbn=%02x, esi=%06x, pid=%08x\n", i, sbn, esi, pid);
		test_assert(sbn == rq_pid2sbn(pid), "%i: sbn:%02x == %02x", i, sbn, rq_pid2sbn(pid));
		test_assert(esi == rq_pid2esi(pid), "%i: esi:%06x == %06x", i, esi, rq_pid2esi(pid));
		pid = 0;
		pid = rq_pidsetsbn(pid, sbn);
		pid = rq_pidsetesi(pid, esi);
		test_assert(sbn == rq_pid2sbn(pid), "%i: sbn:%02x == %02x", i, sbn, rq_pid2sbn(pid));
		test_assert(esi == rq_pid2esi(pid), "%i: esi:%06x == %06x", i, esi, rq_pid2esi(pid));
	}

	return fails;
}
