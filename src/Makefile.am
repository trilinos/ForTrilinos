# @HEADER
# ************************************************************************
# 
#                     
#                 Copyright () Sandia Corporation
# 
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
# 
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#  
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#  
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
# 
# ************************************************************************
# @HEADER

# The following line helps the test harness recover from build errors.

all-local:
	@echo ""
	@echo "Trilinos package ForTrilinos subdirectory src built successfully."
	@echo ""

# Automake does not know about the .f03 extension

.f03.o:
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c -o $@ $<

.f03.obj:
	$(FC) $(AM_FCFLAGS) $(FCFLAGS) -c -o $@ `$(CYGPATH_W) '$<'`

# ------------------------------------------------------------------------
# For each category, create two variables - NAME and NAME_H. The
# second is the list of headers to be installed, i.e., any header that
# might someday be needed by some other code outside New_Package. The first is
# the list of all source and any other header files.
# ------------------------------------------------------------------------

CORE_H =

CORE = \
  forepetra.f03 \
  Epetra_Map_Module.f03 \
  Epetra_Vector_Module.f03 \
  Epetra_SerialComm_Module.f03 \
  Epetra_MpiComm_Module.f03

UTIL_H =

UTIL = 

#AM_CPPFLAGS = -I$(top_srcdir)/../epetra/src -I$(top_builddir)/../epetra/src

# ------------------------------------------------------------------------
# ForTrilinos library specifications
# ------------------------------------------------------------------------
FORTRILINOS_LIB = libfortrilinos.a

FORTRILINOS_H = \
	$(CORE_H) \
	$(UTIL_H)

libfortrilinos_a_SOURCES = \
	$(CORE) \
	$(UTIL)

#EXTRA_libfortrilinos_a_SOURCES =

#We need the following to support separate build directories
# ------------------------------------------------------------------------
# For using a special archiver
# ------------------------------------------------------------------------

# SPECIAL NOTE: Why is the "else" there?  For some reason, automake
# fails to define <lib>_AR of USE_ALTERNATE_AR is false, so we have to
# define it manually.

if USE_ALTERNATE_AR

libfortrilinos_a_AR = $(ALTERNATE_AR)

else

libfortrilinos_a_AR = $(AR) cru

endif

# ------------------------------------------------------------------------
# Some C++ compilers create extra .o-files for templates. We need to
# be sure to include these, and this is the hack to do it.
# ------------------------------------------------------------------------

libfortrilinos_a_LIBADD = $(XTRALDADD)

# ------------------------------------------------------------------------
# List of all libraries to install in $(libexecdir)
# ------------------------------------------------------------------------

lib_LIBRARIES = $(FORTRILINOS_LIB)

# ------------------------------------------------------------------------
# List of all headers to install in $(includedir)
# ------------------------------------------------------------------------

include_HEADERS = $(FORTRILINOS_H) 

# ------------------------------------------------------------------------
# Special stuff to install in our special $(execincludedir)
# ------------------------------------------------------------------------

# SPECIAL NOTE: ForTrilinos_config.h is a machine-dependent file, so we need
# to install it in the machine-dependent directory. However, that is
# not a default installation directory, so we had to create it
# special.

# All Trilinos headers are now installed in the same directory
execincludedir = $(exec_prefix)/include
#This file does not exist at this time
#nodist_execinclude_HEADERS = ForTrilinos_config.h

# ------------------------------------------------------------------------
# Files to be deleted by 'make maintainer-clean'
# ------------------------------------------------------------------------

MAINTAINERCLEANFILES = Makefile.in 

