# Process this file with autoconf to produce a configure script.
AC_INIT(rgraph, [2.2.1], [roger.guimera@gmail.com])
AC_REVISION($Revision$)
AM_INIT_AUTOMAKE
LT_INIT
AC_DISABLE_STATIC
AC_PROG_LIBTOOL

# Create the config.h.in
AC_CONFIG_SRCDIR([config.h.in])
AC_CONFIG_HEADER(config.h)
AC_CONFIG_MACRO_DIR([m4])

# Checks for programs.
AC_PROG_CC
gl_EARLY
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET

AC_PROG_MAKE_SET

# Checks for libraries.
AC_CHECK_LIB([m], [exp])
AC_CHECK_LIB([gslcblas], [cblas_zgemm], , AC_MSG_ERROR([GSL CBLAS library not found]))
AC_CHECK_LIB([gsl], [gsl_vector_set], , AC_MSG_ERROR([GSL library not found]))

# For gnulib
gl_INIT

# Checks for header files.
AC_HEADER_STDC
AC_HEADER_STDBOOL
AC_CHECK_HEADERS([stdlib.h string.h])
AC_CHECK_HEADERS([gsl/gsl_vector.h gsl/gsl_rng.h], , AC_MSG_ERROR([Necessary GSL headers not found]))

# Checks for typedefs, structures, and compiler characteristics.
AC_TYPE_SIZE_T
AC_C_CONST

# Checks for library functions.
AC_FUNC_MALLOC
AC_CHECK_FUNCS([floor pow sqrt])

AC_CONFIG_FILES(lib/Makefile
		src/Makefile
		test/Makefile
		netcarto/Makefile
		utils/Makefile
		reliability/Makefile
		recommender/Makefile
		conflict/Makefile
		multiblock/Makefile
		only_deg/Makefile
		Makefile)

AC_OUTPUT
