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
#ifdef HAVE_LIBSODIUM
# include <sodium.h>
#endif

#define _TESTING 1

extern int fails;

void fail_msg(char *msg, ...);
void test_assert(int condition, char *msg, ...);
void test_assert_s(int condition);
void test_sleep(time_t tv_sec, long tv_nsec);
void test_strcmp(char *str1, char *str2, char *msg, ...);
void test_strncmp(char *str1, char *str2, size_t len, char *msg, ...);
void test_expect(char *expected, char *got);
void test_expectn(char *expected, char *got, size_t len);
void test_expectiov(struct iovec *expected, struct iovec *got);
void test_log(char *msg, ...);
void test_rusage();
void test_name(char *str, ...);
int test_skip(char *str, ...);
#ifdef HAVE_LIBSODIUM
#define test_randombytes randombytes_buf
#else
void test_randombytes(void *buf, size_t len);
#endif
