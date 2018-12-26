.. _install_fortrilinos:

Installation
============

This section provide guidelines for installing ForTrilinos and its TPLs.

Install third-party libraries
-----------------------------

The following third party libraries (TPLs) are used by ForTrilinos:

+------------------------+------------+---------+
| Packages               | Dependency | Version |
+========================+============+=========+
| BLAS/LAPACK            | Required   | N/A     |
+------------------------+------------+---------+
| MPI                    | Required   | N/A     |
+------------------------+------------+---------+

.. warning::

    Currently, the ForTrilinos library can be built using any Fortran compiler.
    However, only gfortran compiler is supported when building tests.

Building ForTrilinos
--------------------

ForTrilinos is configured and built using `TriBITS <https://tribits.org>`_. ForTrilinos builds
within Trilinos effectively as an internal package. The following steps are
required to build and install ForTrilinos:

1. Download Trilinos version

  ForTrilinos strives to be able to work with the most recent ``develop``
  version of Trilinos. However, being a separate project, sometimes upstream
  changes lead to incompatibility. However, the most recent version of
  ``develop`` that ForTrilinos works with can be found in
  ``docker/trilinos_version``.

  To download Trilinos, users should download Trilinos repository like this:

  .. code::

      $ git clone https://github.com/trilinos/Trilinos.git $TRILINOS_DIR
      $ cd $TRILINOS_DIR
      $ git checkout develop

  Here, ``$TRILINOS_DIR`` is the name you want give to the repository.

2. Download and link ForTrilinos into the Trilinos packages directory

  .. code::

      $ git clone https://github.com/trilinos/ForTrilinos.git $FORTRILINOS_DIR
      $ ln -s $FORTRILINOS_DIR $TRILINOS_DIR/ForTrilinos

3. Create a build directory and in it a CMake configuration script.

  .. code::

      $ cd $TRILINOS_DIR
      $ mkdir build
      $ cd build

  The following CMake script defines a near minimal CMake configuration to build
  Trilinos with ForTrilinos enabled:

  .. code-block:: bash

      #!/usr/bin/env bash

      EXTRA_ARGS=$@

      ARGS=(
          -D CMAKE_BUILD_TYPE=Debug

          -D BUILD_SHARED_LIBS=ON

          ### COMPILERS AND FLAGS ###
          -D Trilinos_ENABLE_Fortran=ON
          -D CMAKE_CXX_FLAGS="-Wall -Wpedantic"

          ### TPLs ###
          -D TPL_ENABLE_MPI=ON
          -D TPL_ENABLE_BLAS=ON
          -D TPL_ENABLE_LAPACK=ON

          ### ETI ###
          -D Trilinos_ENABLE_EXPLICIT_INSTANTIATION=ON
            -D Tpetra_INST_INT_INT=OFF
            -D Tpetra_INST_INT_LONG=OFF
            -D Tpetra_INST_INT_LONG_LONG=ON

          ### PACKAGES CONFIGURATION ###
          -D Trilinos_ENABLE_ALL_PACKAGES=OFF
          -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF

          -D Trilinos_ENABLE_TESTS=OFF
          -D Trilinos_ENABLE_EXAMPLES=OFF

          -D Trilinos_EXTRA_REPOSITORIES="ForTrilinos"

          -D Trilinos_ENABLE_Amesos2=ON
          -D Trilinos_ENABLE_Anasazi=ON
          -D Trilinos_ENABLE_Belos=ON
          -D Trilinos_ENABLE_Epetra=OFF
          -D Trilinos_ENABLE_Ifpack2=ON
          -D Trilinos_ENABLE_MueLu=ON
          -D Trilinos_ENABLE_Stratimikos=ON
          -D Trilinos_ENABLE_Tpetra=ON
          -D Trilinos_ENABLE_Thyra=ON

          ### FORTRILINOS ###
          -D Trilinos_ENABLE_ForTrilinos=ON
              -D ForTrilinos_ENABLE_EXAMPLES=ON
              -D ForTrilinos_ENABLE_TESTS=ON
          )
      cmake "${ARGS[@]}" $EXTRA_ARGS $TRILINOS_DIR

.. _patches:

4. Run the configuration script from your build directory.  Here the CMake
   configure script is assumed to be named ``do-configure``

  .. code::

      $ cd $TRILINOS_DIR/build
      $ ./do-configure

  More install scripts can be found in ``scripts/`` directory in the ForTrilinos
  source tree.

Build this documentation
------------------------

(Re)configure with ``-D ForTrlinos_ENABLE_ReadTheDocs=ON`` and run:

.. code::

    $ make docs

Open ``index.html`` in ``$TRILINOS_DIR/packages/ForTrilinos/docs/html``.

.. note::

   Building the documentation requires the Sphinx html theme
   ``sphinx_rtd_theme`` which does not come installed by default on some
   installations of Sphinx.  Be sure to install ``sphinx_rtd_theme`` (via
   ``pip``, ``conda``, etc.) before building the documentation or build
   errors will occur.
