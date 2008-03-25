# $LastChangedDate$
# $Revision$

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = netcompare betweenness

netcompare_SOURCES = main_netcompare.c
netcompare_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

betweenness_SOURCES = main_betweenness.c
betweenness_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la
