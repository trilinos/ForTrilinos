
EXEEXT = .exe

#
# Actual programs
#

# Automake does not know about the .f03 extension

.f03.o:
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c -o $@ $<

.f03.obj:
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c -o $@ `$(CYGPATH_W) '$<'`

noinst_PROGRAMS = simpleEpetraExteranCInterface simpleEpetraFortranInterface

simpleEpetraExteranCInterface_SOURCES = \
  simpleEpetraExteranCInterface.c \
  FEpetra_Map.cpp \
  FEpetra_Vector.cpp
simpleEpetraExteranCInterface_DEPENDENCIES = $(common_dependencies)
simpleEpetraExteranCInterface_LDADD = $(common_ldadd)

simpleEpetraFortranInterface_SOURCES = \
  FEpetra_Map.cpp \
  FEpetra_Vector.cpp \
  simpleEpetraFortranInterface.f03
simpleEpetraFortranInterface_DEPENDENCIES = $(common_dependencies)
simpleEpetraFortranInterface_LDADD = $(common_ldadd)

forepetraext.mod: ../src/forepetraext.mod
	ln -sf ../src/forepetraext.mod .

#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/../epetra/src/libepetra.a  forepetraext.mod

include $(top_builddir)/../epetraext/Makefile.export.epetraext

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(EPETRAEXT_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(EPETRAEXT_LIBS))
else
EXPORT_INCLUDES = $(EPETRAEXT_INCLUDES)
EXPORT_LIBS = $(EPETRAEXT_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# This is already added part of EPETRAEXT_LIBS and therefore automake does not need to add this again!
#LIBS =


EXTRA_DIST = definition
