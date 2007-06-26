TESTS = graph1

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src

check_PROGRAMS = graph1
graph1_SOURCES = main_graph1.c
graph1_LDADD = $(top_srcdir)/lib/libgnu.la $(top_srcdir)/src/librgraph.a
