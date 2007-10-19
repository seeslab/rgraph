# $LastChangedDate$
# $Revision$

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src

bin_PROGRAMS = netcarto
netcarto_SOURCES = main_netcarto.c
netcarto_LDADD = $(top_srcdir)/src/librgraph.la $(top_srcdir)/lib/libgnu.la

CLEANFILES = modules.clu modules.dat network.net node_prop.dat \
	randomized_mod.dat roles.clu roles.dat
EXTRA_DIST = testnet.dat
