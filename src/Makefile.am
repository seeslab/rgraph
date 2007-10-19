INCLUDES		=  -I$(top_builddir) -I$(top_srcdir) \
			-I$(top_srcdir)/lib

pkglib_LTLIBRARIES 	= librgraph.la
librgraph_la_LIBADD	= $(top_srcdir)/lib/libgnu.la
librgraph_la_SOURCES 	= tools.c graph.c models.c modules.c bipartite.c

pkginclude_HEADERS 	= tools.h graph.h models.h modules.h bipartite.h
nodist_pkginclude_HEADERS 	= $(top_srcdir)/lib/search.h
