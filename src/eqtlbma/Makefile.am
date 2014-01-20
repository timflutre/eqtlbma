SUBDIRS = utils irls quantgen tabix

lib_LTLIBRARIES = libeqtlbma.la

libeqtlbma_la_SOURCES = \
	utils/utils_io.hpp \
	utils/utils_math.hpp \
	irls/IRLS.h \
	irls/LogLink.h \
	irls/LinkFunc.h \
	quantgen/covariates.hpp \
	quantgen/gene.hpp \
	quantgen/gene_snp_pair.hpp \
	quantgen/grid.hpp \
	quantgen/MVLR.hpp \
	quantgen/samples.hpp \
	quantgen/snp.hpp \
	tabix/bam_endian.h \
	tabix/bgzf.h \
	tabix/khash.h \
	tabix/knetfile.h \
	tabix/kseq.h \
	tabix/ksort.h \
	tabix/kstring.h \
	tabix/tabix.h

# Dummy C++ source to cause C++ linking.
nodist_EXTRA_libeqtlbma_la_SOURCES = dummy.cxx

libeqtlbma_la_LIBADD = \
	utils/libutils.la \
	irls/libirls.la \
	quantgen/libquantgen.la \
	tabix/libtabix.la
