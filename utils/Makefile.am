# $LastChangedDate$
# $Revision$

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = netcompare netprop betweenness countlinks missing netrandomize

netcompare_SOURCES = main_netcompare.c
netcompare_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

netprop_SOURCES = main_netprop.c
netprop_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

betweenness_SOURCES = main_betweenness.c
betweenness_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

countlinks_SOURCES = main_countlinks.c
countlinks_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

missing_SOURCES = main_missing.c
missing_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

netrandomize_SOURCES = main_netrandomize.c
netrandomize_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la
