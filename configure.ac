AC_INIT([eqtlbma], [1.1], [tflutre@uchicago.edu])
AC_CONFIG_SRCDIR([src/eqtlbma.cpp]) # safety check
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror])
AC_PROG_CXX
AC_CONFIG_FILES([Makefile src/Makefile doc/Makefile tests/Makefile scripts/Makefile])
AC_CHECK_LIB(z, deflatePending, [], [AC_MSG_ERROR([unable to find the zlib library with version >= 1.2.6])])
AC_CHECK_LIB(gslcblas, cblas_sgemm, [], [AC_MSG_ERROR([unable to find the gslcblas library with the cblas_sgemm() function])], [-lgsl]cblas)
AC_CHECK_LIB(gsl, gsl_stats_tss, [], [AC_MSG_ERROR([unable to find the gsl library with the gsl_stats_tss() function])], [-lgsl])
AC_OUTPUT
