/* SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only */
/* Copyright (c) 2022 Brett Sheffield <bacs@librecast.net> */

/* gather - example program to read symbols from separate files created by
 * scatter and reconstitute into original file */

#include <lcrq.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <arpa/inet.h>
#include <fcntl.h>
#include <dirent.h>
#include <libgen.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>

#define T 1024 /* symbol size */

static int decode_symbols(const char *file, size_t F, mode_t mode, uint8_t *enc, uint32_t ESI[], int nesi)
{
	rq_t *rq;
	uint8_t *dec;
	ssize_t byt;
	int rc = 0, fd;

	/* decode symbols */
	rq = rq_init(F, T);
	dec = calloc(rq_KP(rq), T);
	rc = rq_decode_block_rfc(rq, dec, enc, ESI, nesi);
	if (rc == 0) {
		/* write decoded data to file */
		fd = open(file, O_CREAT | O_EXCL | O_WRONLY, mode);
		if (fd != -1) {
			for (; F; F -= byt) {
				byt = write(fd, dec, F);
				if (byt == -1) break;
			}
			close(fd);
		}
	}
	free(dec);

	/* free RaptorQ context */
	rq_free(rq);

	return (F == 0) ? 0 : -1;
}

static int filter_sym(const struct dirent *dent_arthurdent, const char *base,
		size_t *F, rq_pid_t *pid, struct stat *sb)
{
	const char xtra[] = "-0.0.sym";
	size_t dlen, blen;

	/* check file is a regular file and the name is long enough */
	if (stat(dent_arthurdent->d_name, sb) == -1) return 0;
	if ((sb->st_mode & S_IFMT) != S_IFREG) return 0;
	dlen = strlen(dent_arthurdent->d_name);
	blen = strlen(base);
	if (dlen <= blen + strlen(xtra)) return 0;

	/* check basename matches and filename ends in *.sym */
	if (strncmp(dent_arthurdent->d_name, base, blen)) return 0;
	if (strncmp(dent_arthurdent->d_name + dlen - 4, ".sym", 4)) return 0;

	/* now extract filesize F and Payload ID */
	const char *ptr = dent_arthurdent->d_name + blen + 1;
	char *endptr = NULL;
	*F = strtoull(ptr, &endptr, 0);
	*pid = strtoul(endptr + 1, NULL, 0);

	return 1;
}

static int usage(const char *progname, const int rc)
{
	fprintf(stderr, "%s will reassemble a file from RaptorQ symbols created with\n"
			"the companion program scatter.\n\n"
			"Symbols of the form filename.size.pid.sym are read from the\n"
			"current directory. These are then decoded and written to the\n"
			"a new file with the original filename in the current directory.\n"
			"Permissions are set using those from the symbols, which match the\n"
			"permissions of the original file\n\n", progname);
	fprintf(stderr, "Usage: %s filename\n", progname);
	return rc;
}

int main(int argc, char *argv[])
{
	const char *prog = basename(argv[0]);
	const char *file = argv[1];
	uint8_t *enc = NULL;
	uint8_t *ptr;
	uint32_t *ESI;
	size_t F = 0;
	rq_pid_t pid;
	struct stat sb = {0};
	struct dirent *dent;
	DIR *dirp;
	ssize_t byt;
	int c = 0, fd, rc = 1, nesi = 0;

	if (argc != 2) return usage(prog, -1);

	/* scan files, count symbols */
	dirp = opendir(".");
	while ((dent = readdir(dirp))) {
		if (filter_sym(dent, file, &F, &pid, &sb)) nesi++;
	}
	closedir(dirp);

	/* allocate buffer and ESI index */
	enc = malloc(nesi * T);
	ESI = calloc(nesi, sizeof(uint32_t));

	/* read symbols into buffer */
	c = 0;
	ptr = enc;
	dirp = opendir(".");
	while (c < nesi && (dent = readdir(dirp))) {
		if (filter_sym(dent, file, &F, &pid, &sb)) {
#if DEBUG
			printf("%s (%zu, %u)\n", dent->d_name, F, pid);
#endif
			if ((fd = open(dent->d_name, O_RDONLY)) == -1) continue;
			byt = read(fd, ptr, T);
			close(fd);
			/* only increment ptr if correct bytes read */
			if (byt == T) {
				ptr += T;
				ESI[c++] = ntohl(pid << 8);
			}
		}
	}
	closedir(dirp);
	printf("read %i symbols", c);

	if (F) rc = decode_symbols(file, F, sb.st_mode, enc, ESI, nesi);
	if (rc == 0) printf(", %zu bytes written to %s\n", F, file);
	else printf("\n");

	free(ESI);
	free(enc);

	return rc;
}
