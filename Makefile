# SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
# Copyright (c) 2022 Brett Sheffield <bacs@librecast.net>

export VERSION := 0.0.0
export ABIVERS := 0.0
PREFIX ?= /usr/local
export PREFIX
export LIBNAME := liblcrq
LIBDIR := $(PREFIX)/lib
LIBFILE := lib${LIBNAME}.so
INCLUDEDIR := $(PREFIX)/include

all: src

install: all doc
	cd src && $(MAKE) $@

uninstall:
	cd src && $(MAKE) $@

.PHONY: clean realclean src test sparse

src:
	$(MAKE) -C $@

fixme:
	grep -n FIXME src/*.{c,h} test/*.{c,h}

todo:
	grep -n TODO src/*.{c,h} test/*.{c,h}

speedtest:
	cd test && $(MAKE) $@

clean realclean:
	cd src && $(MAKE) $@
	cd test && $(MAKE) $@

sparse: clean
	CC=cgcc $(MAKE) src

clang: clean
	CC=clang $(MAKE) src

clangtest: clean
	CC=clang $(MAKE) test

gcc: clean all

check test sanitize: clean src
	cd test && $(MAKE) $@

%.test %.check %.debug: src
	cd test && $(MAKE) $@
