/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <assert.h>
#include <lcrq.h>

int main(void)
{
	loginit();

	test_name("Lets Get This Party Started...");

	test_rusage();

	return fails;
}
