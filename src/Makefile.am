INCLUDES		=  -I$(top_builddir) -I$(top_srcdir) \
			-I$(top_srcdir)/lib

#pkglib_LTLIBRARIES 	= librgraph.la
#librgraph_la_LIBADD	= $(top_srcdir)/lib/libgnu.la
#librgraph_la_SOURCES 	= tools.c graph.c models.c modules.c bipartite.c missing.c

pkglib_LIBRARIES 	= librgraph.a
librgraph_a_LIBADD	= $(top_srcdir)/lib/libgnu.la
librgraph_a_SOURCES 	= datastruct.c tools.c graph.c models.c modules.c bipartite.c missing.c matrix.c layout.c recommend.c conflict.c genetic.c $(top_srcdir)/lib/search.c

pkginclude_HEADERS 	= datastruct.h tools.h graph.h models.h modules.h bipartite.h missing.h matrix.h layout.h recommend.h conflict.h genetic.h $(top_srcdir)/lib/search.h
