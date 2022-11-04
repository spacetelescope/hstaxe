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
AC_CANONICAL_HOST
PKG_PROG_PKG_CONFIG

# libtools init and update
LT_PREREQ([2.4.2])
LT_INIT
AC_SUBST([LIBTOOL_DEPS])

# RPATH is the path used to look for shared library files.
AC_MSG_CHECKING(RPATH)
if test -z "$RPATH"
then
	case $host_os in
		*linux*) RPATH='-Wl,-rpath,\$$ORIGIN/../lib';;
		*darwin*) RPATH='-Wl,-rpath,${prefix}/lib';;
		*) RPATH='';;
	esac
fi
AC_SUBST([RPATH])
AC_MSG_RESULT($RPATH)
AC_MSG_RESULT($host_os)

LDFLAGS="$LDFLAGS -lm $RPATH"


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

AC_ARG_WITH(wcstools,
    [AS_HELP_STRING([--with-wcstools],
                    [path to wcstools])],
    [use_wcstools=$withval],
    [use_wcstools=no])


AC_ARG_WITH([wcstools-libname],
    [AS_HELP_STRING([--with-wcstools-libname],
                    [set wcstools library name])],
    [WCSTOOLS_LIBNAME=$withval],
    [WCSTOOLS_LIBNAME=wcs])


AC_ARG_WITH(wcstools-libname, [AS_HELP_STRING([--with-wcstools-libname], [path to wcstools])])
AS_IF([test "x$with_wcstools_libname" != "xno"], [
       WCSTOOLS_LIBNAME="$withval"
       ])

if test "x$GSL" != "x"; then
    GSL_LIBS=`$GSL_CONFIG --libs`
    GSL_CFLAGS=`$GSL_CONFIG --cflags`
else
    PKG_CHECK_MODULES([GSL], [gsl >= 1])
fi

if test "x$CFITSIO" != "x"; then
    CFITSIO_LIBS="-L$CFITSIO/lib -lcfitsio -lpthread"
    CFITSIO_CFLAGS="-I$CFITSIO/include"
else
    PKG_CHECK_MODULES([CFITSIO], [cfitsio >= 3])
fi

if test "x$use_wcstools" == "xno"; then
    AC_MSG_NOTICE([using wcstools from... pkg-config])
    PKG_CHECK_MODULES([WCSTOOLS], [wcstools >= 3.9],
        [
            if $PKG_CONFIG --libs wcstools | grep 'lwcstools' >/dev/null; then
                WCSTOOLS_LIBNAME=wcstools
            fi
        ])
else
    AC_MSG_NOTICE([using wcstools from... user-defined])
    incdir="$WCSTOOLS/include"
    WCSTOOLS_LIBS="-L$WCSTOOLS/lib -l${WCSTOOLS_LIBNAME} -lm"
    WCSTOOLS_CFLAGS="-I$incdir -I$incdir/${WCSTOOLS_LIBNAME}"
fi

LDFLAGS="$LDFLAGS $WCSTOOLS_LIBS"
CFLAGS="$CFLAGS $WCSTOOLS_CFLAGS"

# Checks for header files.
AC_CHECK_HEADERS([limits.h math.h stddef.h stdlib.h string.h unistd.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_CHECK_HEADER_STDBOOL

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor pow rint sqrt strchr strdup strstr strtol])

cat << EOF

hstaxe has been configured with the following

Host:                ${host}
Compiler:            ${CC}
Preprocessor flags:  ${CPPFLAGS}
Compiler flags:      ${CFLAGS}
Preprocessor flags:  ${CPPFLAGS}
Linker flags:        ${LDFLAGS}

CFITSIO_CFLAGS:   $CFITSIO_CFLAGS
CFITSIO_LIBS:     $CFITSIO_LIBS
GSL_CFLAGS:       $GSL_CFLAGS
GSL_LIBS:         $GSL_LIBS
WCSTOOLS_CFLAGS:  $WCSTOOLS_CFLAGS
WCSTOOLS_LIBNAME: $WCSTOOLS_LIBNAME
WCSTOOLS_LIBS:    $WCSTOOLS_LIBS

EOF

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT

