TESTS = graph1 giant coclas bipart1

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src

check_PROGRAMS = graph1 giant coclas bipart1

graph1_SOURCES = main_graph1.c
graph1_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

giant_SOURCES = main_giant.c
giant_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

coclas_SOURCES = main_coclas.c
coclas_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

bipart1_SOURCES = main_bipart1.c
bipart1_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la
