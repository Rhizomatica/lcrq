/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2020-2023 Brett Sheffield <bacs@librecast.net> */

#include "test.h"
#include <fcntl.h>
#include <semaphore.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <unistd.h>
#ifdef HAVE_SYS_RANDOM_H
# include <sys/random.h>
#endif

#define DEFAULT_TERM_COLS 80
#define TESTID_WIDTH 16

int COLS = DEFAULT_TERM_COLS;
int MSGW;
int test_status = TEST_OK; /* test exit status (not a count of failures) */
int capreqd = 0;
int capfail = 0;
sem_t log_lock;

void vfail_msg(char *msg, va_list argp)
{
	char *b;
	b = malloc(_vscprintf(msg, argp) + 1);
	vsprintf(b, msg, argp);
	printf("\n            %-*s", MSGW, b);
	free(b);
	test_status = TEST_FAIL;
}

void fail_msg(char *msg, ...)
{
	va_list argp;
	va_start(argp, msg);
	vfail_msg(msg, argp);
	va_end(argp);
}

int test_assert(int condition, char *msg, ...)
{
	if (!condition) {
		va_list argp;
		va_start(argp, msg);
		vfail_msg(msg, argp);
		va_end(argp);
	}
	return condition;
}

void test_assert_s(int condition)
{
	test_assert(condition, "test");
}

int test_strcmp(char *str1, char *str2, char *msg, ...)
{
	if (str1 == NULL || str2 == NULL || strcmp(str1, str2)) {
		va_list argp;
		va_start(argp, msg);
		vfail_msg(msg, argp);
		va_end(argp);
		return 0;
	}
	return 1;
}

int test_strncmp(char *str1, char *str2, size_t len, char *msg, ...)
{
	if (str1 == NULL || str2 == NULL || strncmp(str1, str2, len)) {
		va_list argp;
		va_start(argp, msg);
		vfail_msg(msg, argp);
		va_end(argp);
		return 0;
	}
	return 1;
}

int test_expect(char *expected, char *got)
{
	return test_strcmp(expected, got, "expected: '%s', got: '%s'", expected, got);
}

int test_expectn(char *expected, char *got, size_t len)
{
	return test_strncmp(expected, got, len, "expected: '%s', got: '%s'", expected, got);
}

void test_log(char *msg, ...)
{
	char *b;
	va_list argp;
	va_start(argp, msg);
	b = malloc(_vscprintf(msg, argp) + 1);
	vsprintf(b, msg, argp);
	sem_wait(&log_lock);
	fprintf(stderr, "%s\n", b);
	sem_post(&log_lock);
	va_end(argp);
	free(b);
}

static void init_terminal(void)
{
#if HAVE_SYS_IOCTL_H
	if (isatty(fileno(stdout))) {
		/* get terminal size */
		struct winsize w;
		COLS = (ioctl(fileno(stdout), TIOCGWINSZ, &w) != -1) ? w.ws_col : DEFAULT_TERM_COLS;
	}
	else COLS = DEFAULT_TERM_COLS;
#endif
	MSGW = COLS - TESTID_WIDTH;
}

void test_name(char *str, ...)
{
	char *b;
	va_list argp;

	init_terminal();
	sem_init(&log_lock, 0, 1);
	if (capfail) {
		printf("%-*s", MSGW, "----- requires capabilities (skipping) -----");
		exit(test_status);
	}
	else if (!capreqd && geteuid() == 0) {
		printf("%-*s", MSGW, "----- does not require root (skipping) -----");
		exit(test_status);
	}
	va_start(argp, str);
	b = malloc(_vscprintf(str, argp) + 1);
	vsprintf(b, str, argp);
	test_log("  (%s)", b);
	printf("%-*s", MSGW, b);
	va_end(argp);
	free(b);
}

int test_skip(char *str, ...)
{
	char *b;
	va_list argp;

	init_terminal();
	sem_init(&log_lock, 0, 1);
	va_start(argp, str);
	b = malloc(_vscprintf(str, argp) + 1);
	vsprintf(b, str, argp);
	printf("(skipped) %-*s", MSGW - 10, b);
	test_log("  (%s)", b);
	va_end(argp);
	free(b);
	return 0;
}

void test_cap_require(int cap)
{
	(void) cap;
	// TODO check for capabilities on Linux
	if (geteuid()) capfail++;
	capreqd++;
}

void test_require_linux(void)
{
#ifndef __linux__
	init_terminal();
	printf("%-*s", MSGW, "----- linux only (skipping) -----");
	exit(test_status);
#endif
}

void test_rusage()
{
	struct rusage ru = {0};
	if (getrusage(RUSAGE_SELF, &ru)) {
		perror("getrusage");
		return;
	}
	test_log("user  :   %lis.%li\n", ru.ru_utime.tv_sec, ru.ru_utime.tv_usec);
	test_log("system:   %lis.%li\n", ru.ru_stime.tv_sec, ru.ru_stime.tv_usec);
	test_log("maxrss:   %li\n", ru.ru_maxrss);
	test_log("ixrss:    %li\n", ru.ru_ixrss);
	test_log("idrss:    %li\n", ru.ru_idrss);
	test_log("isrss:    %li\n", ru.ru_isrss);
	test_log("minflt:   %li\n", ru.ru_minflt);
	test_log("majflt:   %li\n", ru.ru_majflt);
	test_log("nswap:    %li\n", ru.ru_nswap);
	test_log("inblock:  %li\n", ru.ru_inblock);
	test_log("oublock:  %li\n", ru.ru_oublock);
	test_log("msgsnd:   %li\n", ru.ru_msgsnd);
	test_log("msgrcv:   %li\n", ru.ru_msgrcv);
	test_log("nsignals: %li\n", ru.ru_nsignals);
	test_log("nvcsw:    %li\n", ru.ru_nvcsw);
	test_log("nivcsw:   %li\n", ru.ru_nivcsw);
}

#ifndef HAVE_LIBSODIUM
void test_randombytes(void *buf, size_t len)
{
	ssize_t return_ignored;
#if defined(HAVE_GETRANDOM) && defined(HAVE_SYS_RANDOM_H)
	return_ignored = getrandom(buf, len, 0);
#else
	static int f; /* we'll keep the handle until program exit */
	if (!f) f = open("/dev/urandom", O_RDONLY);
	return_ignored = read(f, buf, len);
#endif
	(void)return_ignored;
}

uint32_t test_randomnumber(const uint32_t upper_bound)
{
	return (uint32_t)random() % upper_bound;
}
#endif
