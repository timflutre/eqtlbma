AC_INIT([eqtlbma], [1.3], [eqtlbma-users@googlegroups.com], [eqtlbma], [https://github.com/timflutre/eqtlbma/wiki])
AC_CONFIG_SRCDIR([src/eqtlbma_bf.cpp]) # safety check
AC_CONFIG_AUX_DIR([build-aux])
AM_INIT_AUTOMAKE([-Wall -Werror 1.13])

AC_PROG_CXX
AM_PROG_AR
LT_INIT

AC_SUBST([AM_CXXFLAGS])
AC_SUBST([AM_LDFLAGS],[-lm])
AC_SUBST([LIBS])

AC_CONFIG_FILES([Makefile src/Makefile src/eqtlbma/Makefile src/eqtlbma/utils/Makefile src/eqtlbma/quantgen/Makefile src/eqtlbma/irls/Makefile src/eqtlbma/tabix/Makefile doc/Makefile tests/Makefile scripts/Makefile])
AC_CHECK_LIB(z, deflatePending, [], [AC_MSG_ERROR([unable to find the zlib library with version >= 1.2.6])])
AC_CHECK_LIB(gslcblas, cblas_sgemm, [], [AC_MSG_ERROR([unable to find the gslcblas library with the cblas_sgemm() function])], [-lgsl]cblas)
AC_CHECK_LIB(gsl, gsl_stats_tss, [], [AC_MSG_ERROR([unable to find the gsl library with the gsl_stats_tss() function])], [-lgsl])
AC_OPENMP

AM_SILENT_RULES([yes])
AC_OUTPUT
