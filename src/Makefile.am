INCLUDES		=  -I$(top_builddir) -I$(top_srcdir) \
			-I$(top_srcdir)/lib

pkglib_LIBRARIES 	= librgraph.a
librgraph_a_LIBADD	= $(top_srcdir)/lib/libgnu.la
librgraph_a_SOURCES 	= tools.c graph.c models.c modules.c

pkginclude_HEADERS 	= tools.h graph.h models.h modules.h
nodist_pkginclude_HEADERS 	= $(top_srcdir)/lib/search.h
