# $LastChangedDate: 2007-10-26 17:05:23 -0500 (Fri, 26 Oct 2007) $
# $Revision: 96 $

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = reconstruct

reconstruct_SOURCES = main_reconstruct.c
reconstruct_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

CLEANFILES = net_reconstructed.dat
