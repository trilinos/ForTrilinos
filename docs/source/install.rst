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
      $ git checkout trilinos-release-12.12-1

  Here, ``$TRILINOS_DIR`` is the name you want give to the repository.

2. Download and link ForTrilinos into the Trilinos packages directory

  .. code::

      $ git clone https://github.com/trilinos/ForTrilinos.git $FORTRILINOS_DIR
      $ ln -s $FORTRILINOS_DIR $TRILINOS_DIR/packages

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

          ### PACKAGES CONFIGURATION ###
          -D Trilinos_ENABLE_ALL_PACKAGES=OFF
          -D Trilinos_ENABLE_ALL_OPTIONAL_PACKAGES=OFF

          -D Trilinos_ENABLE_TESTS=OFF
          -D Trilinos_ENABLE_EXAMPLES=OFF

          -D Trilinos_ENABLE_Amesos2=ON
          -D Trilinos_ENABLE_Anasazi=ON
          -D Anasazi_ENABLE_Epetra=OFF
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

4. Apply one-time patches to Trilinos.  In future updates, patches will be incorporated directly in to Trilinos but, for now, they must be applied manually.  All required patches can be found in ``scripts/patches`` directory in the ForTrilinos source tree. The script ``scripts/patches/apply-patches`` can be used to apply all of the patches at once.

  .. code::

     $ cd $FORTRILINOS_DIR/scripts/patches
     $ ./apply-patches

  .. note::

     Note, the environment variable ``TRILINOS_DIR`` must be defined and point to the 12.12 release version of Trilinos cloned in Step 1 for the patching and building to succeed.


5. Run the configuration script from your build directory.  Here the CMake
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
