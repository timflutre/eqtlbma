TESTS = test_basic.bash \
        test_genes-absent-in-some-subgroups.bash
EXTRA_DIST = functional_tests.R \
             test_basic.bash \
             test_genes-absent-in-some-subgroups.bash
AM_TESTS_ENVIRONMENT = \
	eqtlbma_abspath='$(abs_top_builddir)'/src/eqtlbma; \
	export eqtlbma_abspath; \
	Rscript_abspath='$(abs_top_srcdir)'/tests/functional_tests.R; \
	export Rscript_abspath;
