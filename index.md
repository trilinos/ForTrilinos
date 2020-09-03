---
layout: fortrilinos
---

# [](#header-1)Published documents

- [Automated Fortran–C++ Bindings for Large-Scale Scientific Applications (2020)](https://ieeexplore.ieee.org/abstract/document/8745480/)
- [Documenting automated Fortran-C++ bindings with SWIG (2019)](https://www.osti.gov/biblio/1557490)
- [Existing Fortran interfaces to Trilinos in preparation for exascale ForTrilinos development (2017)](files/ForTrilinos_Existing_Interfaces.pdf)
- [ForTrilinos design document (2017)](files/ForTrilinos_Design_Document.pdf)

# [](#header-2)ForTrilinos installation

ForTrilinos is built as an independent software library with [the Trilinos
software library](https://trilinos.github.io/index.html) as its only
dependency. ForTrilinos can be installed through [the Spack HPC package
manager](https://spack.readthedocs.io/en/latest/) or independently from your
local installation of Trilinos.

To install ForTrilinos version `2.0.0-dev2` through an existing Spack
installation (v0.16 or higher, or the `develop` branch):
```console
$ spack install fortrilinos@2.0.0-dev2 ^trilinos@12.18.1+nox+stratimikos
```

To install manually, you can point to the target Trilinos installation with the
`CMAKE_PREFIX_PATH` or `Trilinos_ROOT` environment variables, or with a
`Trilinos_ROOT` CMake variable:
```console
$ git clone https://github.com/trilinos/ForTrilinos && cd ForTrilinos
$ mkdir build && cd build
$ cmake -DTrilinos_ROOT=/opt/trilinos -DCMAKE_INSTALL_PREFIX=/opt/fortrilinos ..
$ make install
```

# [](#header-3)Modules

ForTrilinos includes several Fortran modules (`forteuchos`, `forbelos`,
`fortpetra`) that are thin layers built on top of Trilinos packages. The
`fortrilinos_hl` package is a high-level set of wrappers that exposes
linear, nonlinear, and eigenvalue solvers.

Different Trilinos packages are required for different levels of functionality:

| Package | Module |
| ------- | ---------- |
| Teuchos | [all] | 
| Belos | `forbelos` |
| Tpetra | `fortpetra` |
| Anasazi | `fortrilinos_hl` |
| NOX | `fortrilinos_hl` |
| Stratimikos | `fortrilinos_hl` |
| Thyra | `fortrilinos_hl` |
| Amesos2 | `fortrilinos_hl` (optional) |
| Ifpack2 | `fortrilinos_hl` (optional) |
| MueLu | `fortrilinos_hl` (optional) |

# [](#header-4)Downstream usage

To use ForTrilinos as part of your CMake app, simply ensure that CMake can find
it (using the standard `CMAKE_PREFIX_PATH` or `ForTrilinos_ROOT` environment
variables, or the `ForTrilinos_ROOT` CMake variable); and add
```cmake
find_package(ForTrilinos)
```

An example application that uses ForTrilinos and MPI-provided Fortran bindings
might look like:
```cmake
cmake_minimum_required(VERSION 3.12)
project(ForTrilinosInstallTest VERSION 0.0.1 LANGUAGES Fortran)

find_package(ForTrilinos)

add_executable(downstream-app downstream-app.F90)
target_link_libraries(downstream-app ForTrilinos::ForTrilinos MPI::MPI_Fortran)
```
and the `downstream-app.F90` app will simply need
```fortran
use forteuchos
use fortpetra
```
