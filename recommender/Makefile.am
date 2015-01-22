# $LastChangedDate: 2007-07-02 11:10:53 -0500 (Mon, 02 Jul 2007) $
# $Revision: 73 $

INCLUDES = -I$(top_builddir) -I$(top_srcdir) \
	-I$(top_srcdir)/lib -I$(top_srcdir)/src
AM_LDFLAGS = -static

bin_PROGRAMS = multi_recommend_2 multi_recommend_k multi_recommend_k_gibbs

multi_recommend_2_SOURCES = main_multi_recommend_2.c
multi_recommend_2_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

multi_recommend_k_SOURCES = main_multi_recommend_k.c
multi_recommend_k_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

multi_recommend_k_gibbs_SOURCES = main_multi_recommend_k_gibbs.c
multi_recommend_k_gibbs_LDADD = $(top_srcdir)/src/librgraph.a $(top_srcdir)/lib/libgnu.la

#CLEANFILES = modules_bipart.dat
#EXTRA_DIST = women.dat
