**************
Infrastructure
**************

This section is intended for system administrators and users who need to
install their own copy of ForTrilinos.

.. _install_fortrilinos:

Installation
============

ForTrilinos is built as an independent software library with `the
Trilinos software library <https://trilinos.github.io/index.html>`_ as
its only dependency. ForTrilinos can be installed through `the Spack HPC
package manager <https://spack.readthedocs.io/en/latest/>`_ or
independently from your local installation of Trilinos.

.. _version:

Version compatibility
---------------------

ForTrilinos wrappers are tightly coupled to the Trilinos API, so upstream
changes require new releases of ForTrilinos. Since the wrapper generation is
dependent on SWIG-Fortran capabilities, changes to SWIG may also affect the
ability to rebuild ForTrilinos wrappers. Note that SWIG is always optional;
the version here simply denotes the version used to generate the included
wrappers.

.. _version_table:

.. table:: Version compatibility table for ForTrilinos.

   ===========  ============== ======================
   ForTrilinos  Trilinos       SWIG
   ===========  ============== ======================
   2.0.0        13             4.0.2+fortran
   2.0.0-dev3   12.18.1        4.0.2+fortran
   2.0.0-dev2   12.18.1        4.0.0+fortran+15e6ed59
   2.0.0-dev1   12.17+8a82b322 4.0.0+fortran+15e6ed59
   1.0          12.8.1         ---
   ===========  ============== ======================

In :ref:`the version table above <version_table>`, the ``+fortran`` suffix for
SWIG indicates `the SWIG-Fortran fork <https://github.com/swig-fortran/swig>`.
``+sha`` refers to a specific Git commit that comes after the given version.

The original implementation of the ForTrilinos was developed prior to 2012.
That code is no longer developed and maintained, and is available using the
``trilinos-release-12-8-1`` tag in the ForTrilinos repository and the
corresponding Trilinos 12.8.1 version.

E4S
---

As of this writing, ForTrilinos is distributed as part of the `E4S Project
<https://e4s-project.github.io/index.html>` and should be available as a
pre-built binary.

Spack
-----

To install ForTrilinos version ``2.0.0`` through an existing Spack
installation (v0.16 or higher, or the ``develop`` branch):

.. code:: console

   $ spack install fortrilinos@2.0.0 ^trilinos+nox+stratimikos

Manual
------

To install manually, you can point to the target Trilinos installation
with the ``CMAKE_PREFIX_PATH`` or ``Trilinos_ROOT`` environment
variables, or with a ``Trilinos_ROOT`` CMake variable:

.. code:: console

   $ git clone https://github.com/trilinos/ForTrilinos && cd ForTrilinos
   $ mkdir build && cd build
   $ cmake -DTrilinos_ROOT=/opt/trilinos -DCMAKE_INSTALL_PREFIX=/opt/fortrilinos ..
   $ make install

The build options are most easily viewed using the ``ccmake`` GUI. All
ForTrilinos options are scoped with a ``ForTrilinos_`` prefix. Some options
(e.g. ``ForTrilinos_USE_MPI``) are derived from the Trilinos installation and
cannot be changed.

Documentation
-------------

To build this documentation, (re)configure with ``-D ForTrilinos_DOCS=ON`` and
run:

.. code:: console

    $ make doc

then open ``$BUILD/doc/html/index.html``.

App infrastructure setup
========================

The ForTrilinos installation is optimized for use with the CMake build
system CMake_. To use ForTrilinos as part of your CMake-based Fortran
app, add

.. code:: cmake

   find_package(ForTrilinos)

to the top level of your ``CMakeLists.txt`` file. To ensure that CMake can find
the ForTrilinos installation, append its install prefix to the standard
``CMAKE_PREFIX_PATH`` or ``ForTrilinos_ROOT`` environment variables, or define
the ``ForTrilinos_ROOT`` CMake variable when configuring your script.

.. code:: console

   $ cmake -DForTrilinos_ROOT=/usr/local/fortrilinos ForTrilinosInstallTest

An example application that uses ForTrilinos and MPI-provided Fortran
bindings might look like:

.. code:: cmake

   cmake_minimum_required(VERSION 3.12)
   project(ForTrilinosInstallTest VERSION 0.0.1 LANGUAGES Fortran)

   find_package(ForTrilinos)

   add_executable(downstream-app downstream-app.F90)
   target_link_libraries(downstream-app
      ForTrilinos::ForTrilinos MPI::MPI_Fortran
   )

where the ``downstream-app.F90`` app will simply need to include the ForTrilinos
modules:

.. code:: fortran

   use forteuchos
   use fortpetra

.. _CMake : https://cmake.org
