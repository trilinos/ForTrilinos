
EXEEXT = .exe

#
# Actual programs
#

# Automake does not know about the .f03 extension

.f03.o:
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c -o $@ $<

.f03.obj:
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c -o $@ `$(CYGPATH_W) '$<'`

#noinst_PROGRAMS = simpleEpetraFortranInterface simpleEpetraOOFortranExample
noinst_PROGRAMS = simpleEpetraOOFortranExample

#simpleEpetraFortranInterface_SOURCES = simpleEpetraFortranInterface.f03 dummy.cpp
#simpleEpetraFortranInterface_DEPENDENCIES = $(common_dependencies)
#simpleEpetraFortranInterface_LDADD = $(common_ldadd)

simpleEpetraOOFortranExample_SOURCES = simpleEpetraOOFortranExample.f03 dummy.cpp
simpleEpetraOOFortranExample_DEPENDENCIES = $(common_dependencies)
simpleEpetraOOFortranExample_LDADD = $(common_ldadd)

module_files : forepetraext.mod epetra_map_module.mod epetra_vector_module.mod

forepetraext.mod: ../src/forepetraext.mod
	ln -sf ../src/forepetraext.mod .

epetra_map_module.mod: ../src/epetra_map_module.mod
	ln -sf ../src/epetra_map_module.mod .

epetra_vector_module.mod: ../src/epetra_vector_module.mod
	ln -sf ../src/epetra_vector_module.mod .

#
# Common stuff
#

common_dependencies = \
  $(top_builddir)/src/libfortrilinos.a forepetraext.mod epetra_map_module.mod epetra_vector_module.mod

include $(top_builddir)/Makefile.export.ForTrilinos

if USING_GNUMAKE
EXPORT_INCLUDES = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_incl_paths.pl $(FORTRILINOS_INCLUDES))
EXPORT_LIBS = $(shell $(PERL_EXE) $(top_srcdir)/config/strip_dup_libs.pl $(FORTRILINOS_LIBS))
else
EXPORT_INCLUDES = $(FORTRILINOS_INCLUDES)
EXPORT_LIBS = $(FORTRILINOS_LIBS)
endif

AM_CPPFLAGS = $(EXPORT_INCLUDES)

common_ldadd = $(EXPORT_LIBS)

# This is already added part of FORTRILINOS_LIBS and therefore automake does not need to add this again!
#LIBS =


EXTRA_DIST = definition
