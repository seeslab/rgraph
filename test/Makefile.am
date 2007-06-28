TESTS = graph1 bipart1

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src

check_PROGRAMS = graph1
graph1_SOURCES = main_graph1.c
graph1_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

check_PROGRAMS = bipart1
bipart1_SOURCES = main_bipart1.c
bipart1_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la
