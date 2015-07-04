# $LastChangedDate$
# $Revision$

INCLUDES = -I$(top_builddir) -I$(top_srcdir) -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = netcarto

netcarto_SOURCES = main_netcarto.c
netcarto_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la
