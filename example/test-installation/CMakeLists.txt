#----------------------------------*-CMake-*----------------------------------#
# Copyright 2020 UT-Battelle, LLC
# License-Filename: LICENSE
# SPDX-License-Identifier: BSD-3-Clause
#-----------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.12)
project(ForTrilinosInstallTest VERSION 0.0.1 LANGUAGES Fortran)

#-----------------------------------------------------------------------------#

find_package(ForTrilinos)

#-----------------------------------------------------------------------------#

add_executable(downstream-app downstream-app.F90)
target_link_libraries(downstream-app ForTrilinos::ForTrilinos MPI::MPI_Fortran)

#-----------------------------------------------------------------------------#

