/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

/* scatter - example program to encode a file and write out symbols to disk */

#include <lcrq.h>
#include <fcntl.h>
#include <libgen.h>
#include <limits.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

#define T 1024 /* symbol size */

/* write symbol to file */
static inline int write_symbol(size_t F, rq_pid_t pid, uint8_t *sym, const char *basename, mode_t mode)
{
	char filename[PATH_MAX];
	ssize_t byt;
	int fd;

	snprintf(filename, sizeof filename, "./%s-%zu.%u.sym", basename, F, pid);
#if DEBUG
	puts(filename);
#endif
	fd = open(filename, O_CREAT | O_EXCL | O_WRONLY, mode);
	if (fd == -1) return -1;
	byt = write(fd, sym, T);
	if (byt != T) return -1;
	close(fd);

	return 0;
}

static int scatter_file(const size_t F, uint8_t *map, const char *basename, mode_t mode)
{
	rq_t *rq;
	uint8_t sym[T];
	uint16_t KP;
	int rc = 0, syms = 0;

	/* initialize RaptorQ context */
	rq = rq_init(F, T);

	/* data precoding (generate intermediate symbols */
	rc = rq_encode_data_rfc(rq, map, F);
	if (rc) goto exit_err;

	/* generate K' + overhead symbols and write to disk */
	KP = rq_KP(rq);
	rq_pid_t pid = 0; /* payload id => SBN (8 bits) + 24 bit ESI in network byte order */
	for (uint16_t i = 0; rc == 0 && i < KP + RQ_OVERHEAD; i++) {
		rq_pkt_gen(rq, &pid, sym, RQ_RAND);
		rc = write_symbol(F, pid, sym, basename, mode);
		if (rc == 0) syms++;
	}
	printf("%i symbols written\n", syms);
exit_err:
	/* free RaptorQ context */
	rq_free(rq);

	return rc;
}

static int usage(const char *progname, const int rc)
{
	fprintf(stderr, "%s will read the specified file, encode it with RaptorQ,\n"
			"writing random repair symbols as separate files in the\n"
			"current directory. These can be reassembled into the\n"
			"original file with the companion program gather.\n\n", progname);
	fprintf(stderr, "Usage: `%s filename\n", progname);
	return rc;
}

int main(int argc, char *argv[])
{
	const char *prog = argv[0];
	const char *file = argv[1];
	char *base;
	uint8_t *map;
	struct stat sb = {0};
	int fd;

	if (argc != 2) return usage(prog, -1);

	/* open input file for reading */
	if (stat(file, &sb) == -1) {
		perror("stat");
		return usage(prog, -1);
	}
	if ((sb.st_mode & S_IFMT) != S_IFREG) {
		fprintf(stderr, "Error: '%s' is not a regular file\n", file);
		return usage(prog, -1);
	}
	fd = open(file, O_RDONLY);
	if (fd == -1) {
		perror("open");
		return 1;
	}

	/* map the file */
	map = mmap(NULL, sb.st_size, PROT_READ, MAP_PRIVATE, fd, 0);
	if (map == MAP_FAILED) {
		perror("mmap");
		close(fd);
		return 1;
	}
	base = strdup(file);
	scatter_file(sb.st_size, map, basename(base), sb.st_mode);
	free(base);

	/* unmap & close file */
	munmap(map, sb.st_size);
	close(fd);

	return 0;
}
