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
dependent on SWIG-Fortran capabilities, changes there will also affect the
ability to rebuild ForTrilinos wrappers. (Note that SWIG is always optional;
the version here simply denotes the version used to generate the included
wrappers.)

.. table:: Version compatibility table for ForTrilinos.

   ===========  ============== ======================
   ForTrilinos  Trilinos       SWIG
   ===========  ============== ======================
   2.0.0        13.0           4.0.2+fortran
   2.0.0-dev3   12.18.1        4.0.2+fortran
   2.0.0-dev2   12.18.1        4.0.0+fortran+15e6ed59
   2.0.0-dev1   12.17+8a82b322 4.0.0+fortran+15e6ed59
   1.0          12.8.1         ---
   ===========  ============== ======================

The original implementation of the ForTrilinos was developed prior to 2012.
That code is no longer developed and maintained, and is available using the
``trilinos-release-12-8-1`` tag and the corresponding Trilinos 12.8.1 version.

Spack
-----

To install ForTrilinos version ``2.0.0-dev2`` through an existing Spack
installation (v0.16 or higher, or the ``develop`` branch):

.. code:: console

   $ spack install fortrilinos@2.0.0-dev2 ^trilinos@12.18.1+nox+stratimikos

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
