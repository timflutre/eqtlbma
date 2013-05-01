SUBDIRS = irls quantgen

lib_LTLIBRARIES = libeqtlbma.la

libeqtlbma_la_SOURCES = irls/IRLS.h \
	irls/LogLink.h \
	irls/LinkFunc.h \
	quantgen/covariates.hpp \
        quantgen/gene.hpp \
        quantgen/gene_snp_pair.hpp \
        quantgen/grid.hpp \
        quantgen/MVLR.hpp \
        quantgen/samples.hpp \
        quantgen/snp.hpp \
        quantgen/utils_io.hpp \
        quantgen/utils_math.hpp

# Dummy C++ source to cause C++ linking.
nodist_EXTRA_libeqtlbma_la_SOURCES = dummy.cxx

libeqtlbma_la_LIBADD = \
	irls/libirls.la \
	quantgen/libquantgen.la 
