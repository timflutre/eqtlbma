#MAINTAINERCLEANFILES = Makefile.in launch_functional_tests.in
TESTS = launch_functional_tests.bash
EXTRA_DIST = launch_functional_tests.bash functional_tests.R
AM_TESTS_ENVIRONMENT = \
	eqtlbma_abspath='$(abs_top_builddir)'/src/eqtlbma; \
	export eqtlbma_abspath; \
	Rscript_abspath='$(abs_top_srcdir)'/tests/functional_tests.R; \
	export Rscript_abspath;
