# NagBody configure.ac file - mgpt code
# Copyright (c) 2006-2018  Mario A. Rodriguez-Meza, Ciudad de Mexico.
#
#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([mgpt], [1.0.0], [marioalberto.rodriguezmeza@gmail.com])

AM_INIT_AUTOMAKE(silent-rules)
AM_SILENT_RULES
# The following trajectory is important to find install-sh for instance.
AC_CONFIG_SRCDIR([src/perturbations/mgpt/main.c])
AC_CONFIG_HEADERS([config.h])

AC_PREFIX_DEFAULT(/usr/local/NagBody)

# Checks for programs.
AC_PROG_CC
AC_PROG_INSTALL
AC_PROG_RANLIB

# check OpenMP support
AC_OPENMP
if test "x$ac_cv_prog_cxx_openmp" != "xunsupported"; then
   CXXFLAGS="$CXXFLAGS $OPENMP_CXXFLAGS"
   AC_DEFINE(HAVE_OPENMP, 1 , Define to 1 if you have OpenMP )
else
   AC_MSG_NOTICE([We don't have OpenMP. Multithreaded operation is disabled])
fi

# Tell automake that this is not a GNU package (and doesn't need NEWS,
# AUTHORS, etc. files)
#AM_INIT_AUTOMAKE([foreign])

# Prevent Makefile from re-running automake etc automatically because
# it always uses the wrong version or doesn't use the necessary flags.
AM_MAINTAINER_MODE

# Checks for libraries.

# Checks for typedefs, structures, and compiler characteristics.

# Checks for command-line options

# Checks for header files.
AC_CHECK_HEADERS([stddef.h sys/param.h sys/timeb.h])

AC_CHECK_FUNCS([pow sqrt strchr strdup])

AC_CHECK_HEADER_STDBOOL
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_TYPE_SIZE_T

# Defining variables to be passed on the compiler via DEFS substitution variable
AC_DEFINE(NOGNU, [1], [Description])
#AC_DEFINE_UNQUOTED(NOGNU)


AC_CONFIG_FILES([Makefile
		doc/Makefile
		doc/perturbations/Makefile
		doc/perturbations/mgpt/Makefile
		src/general_libs/Makefile
		src/general_libs/general/Makefile
		src/general_libs/io/Makefile
		src/general_libs/math/Makefile
		src/perturbations/Makefile
		src/perturbations/mgpt/Makefile
		src/Makefile])

AC_OUTPUT

echo \
"-------------------------------------------------------
${PACKAGE_NAME} Version ${PACKAGE_VERSION}
Prefix: '${prefix}'.
compiler: '${CC} ${CFLAGS} ${CPPFLAGS}'

Package features:
  Async Execution: ${async_exec}

Now type 'make @<:@<target>@:>@'
  where the optional <target> is:
    all			- build all binaries
    install		- install everything

-------------------------------------------------------"


