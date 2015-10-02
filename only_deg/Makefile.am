# $LastChangedDate: 2007-10-26 17:05:23 -0500 (Fri, 26 Oct 2007) $
# $Revision: 96 $

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = only_degeneration_mb only_degeneration_mb_gibbs

only_degeneration_mb_SOURCES = main_only_degeneration_mb.c
only_degeneration_mb_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

only_degeneration_mb_gibbs_SOURCES = main_only_degeneration_mb_gibbs.c
only_degeneration_mb_gibbs_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la
