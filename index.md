---
layout: fortrilinos
---

# [](#header-1)Published documents

- [Existing Fortran interfaces to Trilinos in preparation for exascale ForTrilinos development](https://www.osti.gov/scitech/biblio/1356940)
- [ForTrilinos design document](files/ForTrilinos_Design_Document.pdf)

# [](#header-1)ForTrilinos installation

The following third party libraries (TPLs) are used by ForTrilinos:

| Packages               | Dependency |
|:---------------------- |:---------- |
| BLAS/LAPACK            | Required   |
| MPI                    | Required   |

ForTrilinos is configured and built using [TriBITS](https://tribits.org). ForTrilinos builds
within Trilinos effectively as an internal package. First, link ForTrilinos into the
Trilinos main directory:

```bash
$ cd $TRILINOS_DIR/packages
$ ln -s $ForTrilinos_DIR
```

Create a `do-configure` script such as:

```bash
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

    -D Trilinos_ENABLE_Amesos=ON
    -D Trilinos_ENABLE_AztecOO=ON
    -D Trilinos_ENABLE_Belos=ON
    -D Trilinos_ENABLE_Epetra=ON
    -D Trilinos_ENABLE_EpetraExt=ON
    -D Trilinos_ENABLE_Ifpack=ON
    -D Trilinos_ENABLE_Ifpack2=ON
    -D Trilinos_ENABLE_Stratimikos=ON
    -D Trilinos_ENABLE_Tpetra=ON
    -D Trilinos_ENABLE_Thyra=ON

    ### FORTRILINOS ###
    -D Trilinos_ENABLE_CTrilinos=ON
    -D Trilinos_ENABLE_ForTrilinos=ON
        -D ForTrilinos_ENABLE_EXAMPLES=ON
        -D ForTrilinos_ENABLE_TESTS=ON
    )
cmake "${ARGS[@]}" $EXTRA_ARGS $TRILINOS_DIR
```

and run it from your build directory:

```bash
    $ mkdir build && cd build
    $ ../do-configure
```
