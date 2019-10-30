#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.69])
AC_INIT([aXe], [3.0], [http://hsthelp.stsci.edu])

AXE_MAJOR_VERSION=3
AXE_MID_VERSION=0
AXE_MINOR_VERSION=0

AC_CONFIG_AUX_DIR(config)
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIRS([m4])

AM_INIT_AUTOMAKE(foreign)

# Checks for programs.
AC_PROG_CC
AC_PROG_LIBTOOL
AC_CANONICAL_HOST

AM_PROG_LIBTOOL

# libtools init and update
LT_PREREQ([2.4.6])
LT_INIT()
AC_SUBST([LIBTOOL_DEPS])

# RPATH is the path used to look for shared library files.
AC_MSG_CHECKING(RPATH)
if test -z "$RPATH"
then
	case $host_os in
	      *linux*) RPATH='-Xlinker -rpath \$ORIGIN/../lib';;
        *darwin*) RPATH='-Wl,-rpath,${prefix}/lib';;
	      *)	RPATH='';;
	esac
fi
AC_SUBST([RPATH])
AC_MSG_RESULT($RPATH)
AC_MSG_RESULT($host_os)

LDFLAGS="$LDFLAGS $RPATH"


# Checks for shared libraries.
AC_ARG_WITH(gsl, [AS_HELP_STRING([--with-gsl], [path to gsl])])
AS_IF([test "x$with_gsl" != "xno"], [
       GSL="$withval"
       ])

AC_PATH_PROG([GSL_CONFIG], [gsl-config], [notfound],[$with_gsl:$PATH])	
if test "$GSL_CONFIG" == "notfound"
  then
    AC_MSG_ERROR([gsl-config not found. See ./configure --help to know how to set a proper gsl path])
  fi
CFLAGS="$CFLAGS `$GSL_CONFIG --cflags`"
LDFLAGS_GSL="$LDFLAGS_INI `$GSL_CONFIG --libs`"


AC_CHECK_HEADERS(
    [gsl/gsl_vector.h gsl/gsl_linalg.h gsl/gsl_errno.h gsl/gsl_spline.h],
    [AC_MSG_RESULT([check for gsl include files... yes])],
    [AC_MSG_ERROR([Cannot find some needed GSL include files])]
)

AC_ARG_WITH(cfitsio, [AS_HELP_STRING([--with-cfitsio], [path to cfitsio])])
AS_IF([test "x$with_cfitsio" != "xno"], [
       CFITSIO="$withval"
       ])

AC_ARG_WITH(wcstools, [AS_HELP_STRING([--with-wcstools], [path to wcstools])])
AS_IF([test "x$with_wcstools" != "xno"], [
       if test ! -d "$withval"; then
           AC_MSG_ERROR([--with-wcstools requires a directory prefix])
       fi

       WCSTOOLS="$withval"
       ])

if test "x$GSL" != "x"; then
    GSL_LIBS=`$GSL/bin/gsl-config --libs`
    GSL_CFLAGS=`$GSL/bin/gsl-config --cflags`
else
    PKG_CHECK_MODULES([GSL], [gsl >= 1])
fi

if test "x$CFITSIO" != "x"; then
    CFITSIO_LIBS="-L$CFITSIO/lib -lcfitsio -lpthread"
    CFITSIO_CFLAGS="-I$CFITSIO/include"
else
    PKG_CHECK_MODULES([CFITSIO], [cfitsio >= 3])
fi

if test "x$WCSTOOLS" != "x"; then
    WCSTOOLS_LIBS="-L$WCSTOOLS/lib -lwcs"
    WCSTOOLS_CFLAGS="-I$WCSTOOLS/include"
    LDFLAGS="$LDFLAGS $WCSTOOLS_LIBS"
    CFLAGS="$CFLAGS $WCSTOOLS_CFLAGS"
else
    AC_MSG_ERROR([wcstools does not support pkg-config; use --with-wcstools=/path/to/wcstools])
fi


# Checks for header files.
AC_CHECK_HEADERS([limits.h stddef.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow rint sqrt strchr strdup strstr strtol])

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT()

