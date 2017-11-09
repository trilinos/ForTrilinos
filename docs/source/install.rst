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

1. Download Trilinos 12.12 release

  ForTrilinos release model requires users to work with a specific version of
  Trilinos, as otherwise the interfaces may have changed. Users may download
  Trilinos 12.12.1 source `here
  <https://github.com/trilinos/Trilinos/archive/trilinos-release-12-12-1.tar.gz>`_.
  An alternative is to checkout a specific version of Trilinos repository like
  this:

  .. code::

      $ git clone https://github.com/trilinos/Trilinos.git $TRILINOS_DIR
      $ cd $TRILINOS_DIR
      $ git checkout trilinos-release-12.12

  Here, ``$TRILINOS_DIR`` is the name you want give to the repository.

2. Download and link ForTrilinos into the Trilinos packages directory

  .. code::

      $ git clone https://github.com/trilinos/ForTrilinos.git $FORTRILINOS_DIR
      $ cd $TRILINOS_DIR/packages
      $ ln -s $FORTRILINOS_DIR

3. Create a CMake configuration script

  Here is an example of the ``do-configure`` configuration script:

  .. code-block:: bash

      #!/usr/bin/env/bash

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

          ### PACKAGES CONFIGURATION ###
          -D Trilinos_ENABLE_ALL_PACKAGES=OFF
          -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF

          -D Trilinos_ENABLE_TESTS=OFF
          -D Trilinos_ENABLE_EXAMPLES=OFF

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

4. Run the configuration script

  From your build directory:

  .. code::

      $ mkdir build && cd build
      $ ./do-configure

  More install scripts can be found in ``scripts/`` directory in the ForTrilinos
  source tree.

Build this documentation
------------------------

(Re)configure with ``-D ForTrlinos_ENABLE_ReadTheDocs=ON`` and run:

.. code::

    $ make docs

Open the ``index.html`` in the directory ``ReadTheDocs/docs/html``.
