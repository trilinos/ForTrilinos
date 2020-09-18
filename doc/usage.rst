Usage
=====

To use ForTrilinos as part of your CMake app, simply ensure that CMake
can find it (using the standard ``CMAKE_PREFIX_PATH`` or
``ForTrilinos_ROOT`` environment variables, or the ``ForTrilinos_ROOT``
CMake variable); and add

.. code:: cmake

   find_package(ForTrilinos)

An example application that uses ForTrilinos and MPI-provided Fortran
bindings might look like:

.. code:: cmake

   cmake_minimum_required(VERSION 3.12)
   project(ForTrilinosInstallTest VERSION 0.0.1 LANGUAGES Fortran)

   find_package(ForTrilinos)

   add_executable(downstream-app downstream-app.F90)
   target_link_libraries(downstream-app ForTrilinos::ForTrilinos MPI::MPI_Fortran)

and the ``downstream-app.F90`` app will simply need

.. code:: fortran

   use forteuchos
   use fortpetra
