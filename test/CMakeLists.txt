#---------------------------------*-CMake-*----------------------------------#
# Copyright 2020 UT-Battelle, LLC
# License-Filename: LICENSE
# SPDX-License-Identifier: BSD-3-Clause
#----------------------------------------------------------------------------#

add_library(ForTrilinosTests
  fortest.F90
)
set_target_properties(ForTrilinosTests
  PROPERTIES Fortran_MODULE_DIRECTORY "${ForTrilinos_Fortran_MODULE_DIRECTORY}"
)
target_link_libraries(ForTrilinosTests PUBLIC ForTrilinos)

if(CMAKE_Fortran_COMPILER_ID STREQUAL GNU)
  target_compile_options(ForTrilinosTests INTERFACE
    $<$<COMPILE_LANGUAGE:Fortran>:-ffree-line-length-none>
  )
endif()

#---------------------------------------------------------------------------##

function(fortrilinos_add_test source)
  get_filename_component(_DIRNAME "${source}" DIRECTORY)
  get_filename_component(_BASENAME "${source}" NAME_WE)
  message("_DIRNAME=${_DIRNAME}")
  message("_BASENAME=${_BASENAME}")
  string(REGEX REPLACE "[^a-zA-Z0-9_]+" "_" _TARGET "${_DIRNAME}/${_BASENAME}")

  return()
  add_executable("${_TARGET}" "${source}")
  set(_EXE_NAME "$<TARGET_FILE:${_TARGET}>")
  add_test(NAME "${_BASENAME}" COMMAND ${_test_cmd} ${_test_args})
endfunction()

#---------------------------------------------------------------------------##

fortrilinos_add_test(forteuchos/test_teuchos_array.F90)
return()
fortrilinos_add_test(forteuchos/test_teuchos_comm.F90)
fortrilinos_add_test(forteuchos/test_teuchos_plist.F90)

if(ForTrilinos_USE_HL)
  fortrilinos_add_test(fortrilinos_hl/test_eigen_handle.cpp)
  fortrilinos_add_test(fortrilinos_hl/test_solver_handle.cpp)
  fortrilinos_add_test(fortrilinos_hl/test_eigen_handle_preset.cpp)
  fortrilinos_add_test(fortrilinos_hl/test_nox.cpp)
endif()

if(TARGET forbelos)
  fortrilinos_add_test(forbelos/test_belos_types.F90)
endif()

if(TARGET fortpetra)
  fortrilinos_add_test(fortpetra/test_tpetra_crsgraph.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_crsgraph_helper.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_crsmatrix.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_crsmatrix_helper.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_import_export.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_import_export_helper.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_map.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_multivector.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_multivector_helper.F90)
endif()

if(ForTrilinos_USE_HL)
  fortrilinos_add_test(fortrilinos_hl/)
  fortrilinos_add_test(fortrilinos_hl/test_eigen_handle.F90)
  fortrilinos_add_test(fortrilinos_hl/test_eigen_unit.F90)
  fortrilinos_add_test(fortrilinos_hl/test_futility_reproducer.F90)
  fortrilinos_add_test(fortrilinos_hl/test_solver_handle.F90)
  fortrilinos_add_test(fortrilinos_hl/test_solver_unit.F90)
  fortrilinos_add_test(fortrilinos_hl/test_nox.F90
    fortrilinos_hl/Tpetra_ModelEvaluator_1DFEM.F90)
endif()

#---------------------------------------------------------------------------##
