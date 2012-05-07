TESTS = graph1 giant coclas bipart1 bipart2 tools betweenness lapspec netsum

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src

check_PROGRAMS = graph1 giant coclas bipart1 bipart2 tools betweenness \
		 lapspec netsum

graph1_SOURCES = main_graph1.c
graph1_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

giant_SOURCES = main_giant.c
giant_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

coclas_SOURCES = main_coclas.c
coclas_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

bipart1_SOURCES = main_bipart1.c
bipart1_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

bipart2_SOURCES = main_bipart2.c
bipart2_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

tools_SOURCES = main_tools.c
tools_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

betweenness_SOURCES = main_betweenness.c
betweenness_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

lapspec_SOURCES = main_lapspec.c
lapspec_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

netsum_SOURCES = main_netsum.c
netsum_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

EXTRA_DIST = test.dat testbinet.dat testbinetpart.dat test_betw.dat test_graph1.dat
CLEANFILES = *.net
