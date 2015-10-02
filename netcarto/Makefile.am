INCLUDES = -I$(top_builddir) -I$(top_srcdir) -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = netcarto netcarto-legacy

netcarto_SOURCES = main_netcarto.c
netcarto_LDADD = $(top_srcdir)/src/librgraph.a 

netcarto_legacy_SOURCES = main_netcarto_legacy.c
netcarto_legacy_LDADD =  $(top_srcdir)/src/librgraph.a
