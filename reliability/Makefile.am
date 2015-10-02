# $LastChangedDate: 2007-10-26 17:05:23 -0500 (Fri, 26 Oct 2007) $
# $Revision: 96 $

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = reliability_reconstruct reliability_links reliability_links_k reliability_links_gibbs reliability_links_gibbs_sparse

reliability_reconstruct_SOURCES = main_reliability_reconstruct.c
reliability_reconstruct_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

reliability_links_SOURCES = main_reliability_links.c
reliability_links_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

reliability_links_k_SOURCES = main_reliability_links_k.c
reliability_links_k_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

reliability_links_gibbs_SOURCES = main_reliability_links_gibbs.c
reliability_links_gibbs_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

reliability_links_gibbs_sparse_SOURCES = main_reliability_links_gibbs_sparse.c
reliability_links_gibbs_sparse_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

CLEANFILES = net_reconstructed.dat bogus.dat missing.dat
