AM_CXXFLAGS = -I$(top_srcdir)/src -fopenmp @AM_CXXFLAGS@

lib_LIBRARIES = libquantgen.a

libquantgen_adir = $(includedir)/quantgen

libquantgen_a_HEADERS = covariates.hpp \
	gene.hpp \
	gene_snp_pair.hpp \
	grid.hpp \
	MVLR.hpp \
	samples.hpp \
	snp.hpp \
	utils_io.hpp \
	utils_math.hpp

libquantgen_a_SOURCES = $(libmylib_a_HEADERS) \
	covariates.cpp \
	gene.cpp \
	gene_snp_pair.cpp \
	grid.cpp \
	MVLR.cpp \
	samples.cpp \
	snp.cpp \
	utils_io.cpp \
	utils_math.cpp
