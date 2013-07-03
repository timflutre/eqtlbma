AM_CXXFLAGS = -I$(top_srcdir)/src/eqtlbma -fopenmp @AM_CXXFLAGS@

bin_PROGRAMS = eqtlbma_bf eqtlbma_hm eqtlbma_avg_bfs

eqtlbma_bf_SOURCES = eqtlbma_bf.cpp
eqtlbma_bf_CPPFLAGS = -fopenmp
eqtlbma_bf_LDFLAGS = -fopenmp #-pg
eqtlbma_bf_LDADD = $(top_builddir)/src/eqtlbma/libeqtlbma.la $(AM_LDFLAGS)

eqtlbma_hm_SOURCES = eqtlbma_hm.cpp hm_classes.hpp hm_methods.cpp
eqtlbma_hm_CPPFLAGS = -fopenmp
eqtlbma_hm_LDFLAGS = -fopenmp #-pg
eqtlbma_hm_LDADD = $(top_builddir)/src/eqtlbma/libeqtlbma.la $(AM_LDFLAGS)

eqtlbma_avg_bfs_SOURCES = eqtlbma_avg_bfs.cpp
eqtlbma_avg_bfs_CPPFLAGS = -fopenmp
eqtlbma_avg_bfs_LDFLAGS = -fopenmp #-pg
eqtlbma_avg_bfs_LDADD = $(top_builddir)/src/eqtlbma/libeqtlbma.la $(AM_LDFLAGS)

EXTRA_DIST = simul_flutre_et_al.cpp
