#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.71])
AC_INIT([lcrq], [0.0.0], [bugs@librecast.net])
AC_SUBST(PACKAGE_ABIVERS, 0.0)
AC_COPYRIGHT(Copyright (c) 2022 Brett Sheffield <bacs@librecast.net>
See COPYING for license details.
)
AC_CONFIG_SRCDIR([src/])
AC_CONFIG_HEADERS([src/config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_LN_S

# Checks for libraries.
AC_CHECK_LIB([dl], [dlsym])
#AC_CHECK_LIB([sodium], [sodium_init])
AC_ARG_WITH(sodium, "use libsodium for encryption and hashing",
            [AC_SUBST([HAVE_LIBSODIUM], [])],
            [with_sodium=check])
        LIBSODIUM=
          AS_IF([test "x$with_sodium" != xno],
            [AC_CHECK_LIB([sodium], [sodium_init],
              [AC_SUBST([LIBSODIUM], ["-lsodium"])
               AC_SUBST([HAVE_LIBSODIUM], ["HAVE_LIBSODIUM := 1"])
               AC_DEFINE([HAVE_LIBSODIUM], [1],
                         [Define if you have libsodium])
              ],
              [if test "x$with_sodium" != xcheck; then
                 AC_MSG_FAILURE([--with-sodium was given, but libsodium not found])
               fi
              ], -lsodium)])


# Checks for header files.
AC_CHECK_HEADERS([arpa/inet.h fcntl.h stdint.h sys/param.h sys/socket.h
		  sys/time.h unistd.h],[],AC_MSG_ERROR([required header file missing]))

# Checks for typedefs, structures, and compiler characteristics.
AC_C_BIGENDIAN
AC_C_INLINE
AC_TYPE_SIZE_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_MMAP
AC_FUNC_REALLOC
AC_CHECK_FUNCS([clock_gettime floor memset munmap])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 test/Makefile])
AC_OUTPUT
