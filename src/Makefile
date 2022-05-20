# SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
# Copyright (c) 2022 Brett Sheffield <bacs@librecast.net>

SHELL = /bin/sh
ifndef VERSION
$(error VERSION not set)
endif
ifndef ABIVERS
$(error ABIVERS not set)
endif
ifndef LIBNAME
$(error LIBNAME not set)
endif
.SUFFIXES:
.SUFFIXES: .c .o
CFLAGS += -O3 -Wall -Wextra -Wpedantic -g -fPIC -I. -I../include
HEADERS = ../include/$(LIBNAME).h $(sort $(wildcard ../include/$(LIBNAME)/*.h))
INSTALL ?= install
LDCONFIG ?= ldconfig
INSTALL_DATA := $(INSTALL) -m 644
PREFIX ?= /usr/local
LIBDIR := $(DESTDIR)$(PREFIX)/lib
INCLUDEDIR := $(DESTDIR)$(PREFIX)/include
OBJECTS := lcrq.o

LIBS :=
SOFILES := lib$(LIBNAME).so
SOVFILES := $(foreach f,$(SOFILES),$f.$(VERSION))
ABIFILES := $(foreach f,$(SOFILES),$f.$(ABIVERS))

all: $(SOFILES)

lib$(LIBNAME).so: $(OBJECTS)
	$(CC) -shared -O3 -fPIC -Wl,-soname,lib$(LIBNAME).so.$(ABIVERS) -o $@ $^ $(LIBS)

%.o: %.c %.h $(HEADERS)
	$(CC) $(CFLAGS) -c $<

install: $(SOFILES)
	$(INSTALL) -d $(LIBDIR)
	$(INSTALL) -d $(INCLUDEDIR)
	$(INSTALL_DATA) lib$(LIBNAME).so $(LIBDIR)/lib$(LIBNAME).so.$(VERSION)
	cp -r ../include/* $(INCLUDEDIR)
	$(LDCONFIG) -n $(LIBDIR)
	ln -sf lib$(LIBNAME).so.$(ABIVERS) $(LIBDIR)/lib$(LIBNAME).so

.PHONY: clean realclean uninstall

uninstall:
	rm -f $(foreach f,$(SOFILES),$(DESTDIR)$(LIBDIR)/$f.$(VERSION))
	$(LDCONFIG) -n $(LIBDIR)
	cd $(LIBDIR)/ && rm -f $(SOFILES)

clean:
	$(RM) *.o $(SOFILES) $(ABIFILES)

realclean: clean
