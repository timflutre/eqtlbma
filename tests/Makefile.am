TESTS = test_basic.bash \
	test_genes-absent-in-some-subgroups.bash \
	test_with-covariates.bash \
	test_mvlr.bash \
	test_common-uniq-inds.bash \
	test_hm.bash \
	test_poisson.bash
EXTRA_DIST = functional_tests.R \
	test_basic.bash \
	test_genes-absent-in-some-subgroups.bash \
	test_with-covariates.bash \
	test_mvlr.bash \
	test_common-uniq-inds.bash \
	test_hm.bash \
	test_poisson.bash
AM_TESTS_ENVIRONMENT = \
	bf_abspath='$(abs_top_builddir)'/src/eqtlbma_bf; \
	export bf_abspath; \
	hm_abspath='$(abs_top_builddir)'/src/eqtlbma_hm; \
	export hm_abspath; \
	Rscript_abspath='$(abs_top_srcdir)'/tests/functional_tests.R; \
	export Rscript_abspath;
