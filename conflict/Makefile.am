# $LastChangedDate: 2007-07-02 11:10:53 -0500 (Mon, 02 Jul 2007) $
# $Revision: 73 $

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = conflict

conflict_SOURCES = main_conflict.c
conflict_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

#CLEANFILES = modules_bipart.dat
#EXTRA_DIST = women.dat
