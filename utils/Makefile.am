# $LastChangedDate$
# $Revision$

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = netcompare netprop nodeprop countlinks missing netrandomize \
	modularbipart lapspec netlayout multinetlayout

netcompare_SOURCES = main_netcompare.c
netcompare_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

netprop_SOURCES = main_netprop.c
netprop_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

nodeprop_SOURCES = main_nodeprop.c
nodeprop_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

countlinks_SOURCES = main_countlinks.c
countlinks_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

missing_SOURCES = main_missing.c
missing_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

netrandomize_SOURCES = main_netrandomize.c
netrandomize_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

modularbipart_SOURCES = main_modularbipart.c
modularbipart_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

lapspec_SOURCES = main_lapspec.c
lapspec_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

netlayout_SOURCES = main_netlayout.c
netlayout_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

multinetlayout_SOURCES = main_multinetlayout.c
multinetlayout_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la
