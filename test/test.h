/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2020-2022 Brett Sheffield <bacs@librecast.net> */

#include "../src/config.h"
#include <assert.h>
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <time.h>
#include "log.h"
#include "misc.h"
#include "valgrind.h"
#ifdef __linux__
#include <linux/capability.h>
#else
#define CAP_NET_ADMIN 12
#endif
#ifdef HAVE_LIBSODIUM
# include <sodium.h>
#endif

enum {
	TEST_ERR  = -1,
	TEST_OK   = 0,
	TEST_WARN = 1,
	TEST_FAIL = 2,
	TEST_UNKN = 3
};

extern int test_status;

void fail_msg(char *msg, ...);
int test_assert(int condition, char *msg, ...);
void test_assert_s(int condition);
int test_strcmp(char *str1, char *str2, char *msg, ...);
int test_strncmp(char *str1, char *str2, size_t len, char *msg, ...);
int test_expect(char *expected, char *got);
int test_expectn(char *expected, char *got, size_t len);
void test_log(char *msg, ...);
void test_rusage();
void test_name(char *str, ...);
int test_skip(char *str, ...);
void test_cap_require(int cap);
void test_require_linux(void);
#ifdef HAVE_LIBSODIUM
# define test_randombytes randombytes_buf
# define test_randomnumber randombytes_uniform
#else
void test_randombytes(void *buf, size_t len);
uint32_t test_randomnumber(const uint32_t upper_bound);
#endif
