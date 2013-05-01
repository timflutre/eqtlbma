AM_CXXFLAGS = -I$(top_srcdir)/src/eqtlbma -fopenmp @AM_CXXFLAGS@

bin_PROGRAMS = bf hm

bf_SOURCES = bf.cpp
bf_CPPFLAGS = -fopenmp
bf_LDFLAGS = -fopenmp #-pg
bf_LDADD = $(top_builddir)/src/eqtlbma/libeqtlbma.la $(AM_LDFLAGS)

hm_SOURCES = hm.cpp hm_classes.h hm_methods.cpp
hm_CPPFLAGS = 
hm_LDFLAGS = #-pg
hm_LDADD = $(top_builddir)/src/eqtlbma/libeqtlbma.la $(AM_LDFLAGS)
