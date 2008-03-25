INCLUDES		=  -I$(top_builddir) -I$(top_srcdir) \
			-I$(top_srcdir)/lib

#pkglib_LTLIBRARIES 	= librgraph.la
#librgraph_la_LIBADD	= $(top_srcdir)/lib/libgnu.la
#librgraph_la_SOURCES 	= tools.c graph.c models.c modules.c bipartite.c

pkglib_LIBRARIES 	= librgraph.a
librgraph_a_LIBADD	= $(top_srcdir)/lib/libgnu.la
librgraph_a_SOURCES 	= datastruct.c tools.c graph.c models.c modules.c bipartite.c

pkginclude_HEADERS 	= datastruct.h tools.h graph.h models.h modules.h bipartite.h
nodist_pkginclude_HEADERS 	= $(top_srcdir)/lib/search.h
