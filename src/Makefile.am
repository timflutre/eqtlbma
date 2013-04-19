AM_CXXFLAGS = -I$(top_srcdir)/src @AM_CXXFLAGS@

bin_PROGRAMS = eqtlbma hm

eqtlbma_SOURCES = eqtlbma.cpp
eqtlbma_CPPFLAGS = -fopenmp
eqtlbma_LDFLAGS = -fopenmp #-pg
eqtlbma_LDADD = $(top_builddir)/src/quantgen/libquantgen.a $(AM_LDFLAGS)

hm_SOURCES = hm.cpp hm_classes.h hm_methods.cpp
hm_CPPFLAGS = 
hm_LDFLAGS = #-pg
hm_LDADD = $(top_builddir)/src/quantgen/libquantgen.a $(AM_LDFLAGS)
