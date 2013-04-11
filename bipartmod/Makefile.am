# $LastChangedDate$
# $Revision$

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = bipartmod bipartmod_w

bipartmod_SOURCES = main_bipartmod.c
bipartmod_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

bipartmod_w_SOURCES = main_bipartmod_w.c
bipartmod_w_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

CLEANFILES = modules_bipart.dat
EXTRA_DIST = women.dat
