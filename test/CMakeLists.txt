#---------------------------------*-CMake-*----------------------------------#
# Copyright 2020 UT-Battelle, LLC
# License-Filename: LICENSE
# SPDX-License-Identifier: BSD-3-Clause
#----------------------------------------------------------------------------#

#----------------------------------------------------------------------------#
# Test library
#----------------------------------------------------------------------------#

set(_SOURCES fortest.F90)
set(_TEST_DEPLIBS)
if(TARGET fortpetra)
  list(APPEND _SOURCES
    fortpetra/test_tpetra_crsgraph_helper.F90
    fortpetra/test_tpetra_crsmatrix_helper.F90
    fortpetra/test_tpetra_import_export_helper.F90
    fortpetra/test_tpetra_multivector_helper.F90
  )
endif()
if(ForTrilinos_USE_HL)
  list(APPEND _SOURCES fortrilinos_hl/Tpetra_ModelEvaluator_1DFEM.F90)
endif()

add_library(ForTrilinosTests ${_SOURCES})
target_link_libraries(ForTrilinosTests PUBLIC ForTrilinos)
if(ForTrilinos_USE_MPI)
  target_link_libraries(ForTrilinosTests PUBLIC MPI::MPI_Fortran)
endif()
target_include_directories(ForTrilinosTests
  PUBLIC "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>"
)

#----------------------------------------------------------------------------#
# Test files
#----------------------------------------------------------------------------#

file(COPY
  fortrilinos_hl/data/LHS_matrix.mat
  fortrilinos_hl/data/RHS_matrix.mat
  fortrilinos_hl/data/davidson.xml
  fortrilinos_hl/data/nox_params.xml
  fortrilinos_hl/data/stratimikos.xml
  DESTINATION "${CMAKE_CURRENT_BINARY_DIR}"
)

#----------------------------------------------------------------------------#
# Test macros
#----------------------------------------------------------------------------#

set(ForTrilinos_TEST_DEFAULT_NUMPROCS 1 4)

function(fortrilinos_add_test source)
  cmake_parse_arguments(_FTAT
    "DISABLE"
    ""
    "NUMPROCS;SOURCES;LINK_LIBRARIES"
    ${ARGN}
  )

  # Create executable
  get_filename_component(_DIRNAME "${source}" DIRECTORY)
  get_filename_component(_BASENAME "${source}" NAME_WE)
  if(source MATCHES ".cpp$")
    set(_BASENAME "${_BASENAME}/cpp")
  endif()
  string(REGEX REPLACE "[^a-zA-Z0-9_]+" "_" _TARGET "${_DIRNAME}/${_BASENAME}")
  add_executable("${_TARGET}" "${source}" ${_FTAT_SOURCES})
  target_link_libraries("${_TARGET}" ForTrilinosTests ${_FTAT_LINK_LIBRARIES})
  set(_EXE_NAME "$<TARGET_FILE:${_TARGET}>")

  # Set number of processors
  if(NOT _FTAT_NP)
    set(_FTAT_NP ${ForTrilinos_TEST_DEFAULT_NUMPROCS})
  endif()

  # Add tests
  foreach(_np IN LISTS _FTAT_NP)
    if(_np GREATER ForTrilinos_MAX_NUMPROCS)
      continue()
    endif()
    if(ForTrilinos_USE_MPI)
      set(_cmd "${MPIEXEC_EXECUTABLE}" ${MPIEXEC_NUMPROC_FLAG} "${_np}"
        ${MPIEXEC_PREFLAGS} "${_EXE_NAME}" ${MPIEXEC_POSTFLAGS} ${_test_args})
    else()
      set(_cmd "${_EXE_NAME}" ${_test_args})
    endif()

    set(_props)
    if(_np GREATER 1)
      set(_suffix "/${_np}")
      list(APPEND _props PROCESSORS "${_np}")
    else()
      set(_suffix)
    endif()
    if(_FTAT_DISABLE)
      list(APPEND _props DISABLED True)
    endif()

    set(_name "${_DIRNAME}/${_BASENAME}${_suffix}")
    add_test(NAME "${_name}" COMMAND ${_cmd})
    if(_props)
      set_tests_properties("${_name}" PROPERTIES ${_props})
    endif()
  endforeach()
endfunction()

#----------------------------------------------------------------------------#
# Tests
#----------------------------------------------------------------------------#

fortrilinos_add_test(forteuchos/test_teuchos_array.F90)
fortrilinos_add_test(forteuchos/test_teuchos_comm.F90)
fortrilinos_add_test(forteuchos/test_teuchos_plist.F90)

if(TARGET forbelos)
  fortrilinos_add_test(forbelos/test_belos_types.F90)
endif()

if(TARGET fortpetra)
  fortrilinos_add_test(fortpetra/test_tpetra_crsgraph.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_crsmatrix.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_import_export.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_map.F90)
  fortrilinos_add_test(fortpetra/test_tpetra_multivector.F90)
endif()

if(ForTrilinos_USE_HL)
  if(NOT Trilinos_Anasazi_FOUND OR NOT Trilinos_Amesos2_FOUND)
    # Throw from muelu:
    # "Unable to construct any direct solver.Plase enable (TPETRA and AMESOS2)
    # or (EPETRA and AMESOS)"
    set(_needs_anasazi DISABLE)
  endif()

  # The C++ tests indirectly use Trilinos headers and thus require directly
  # linking to Trilinos
  set(_link LINK_LIBRARIES
    Trilinos::Anasazi
    Trilinos::Ifpack2
    Trilinos::NOX
    Trilinos::Stratimikos
    Trilinos::Thyra
    Trilinos::ThyraTpetraAdapters
  )

  fortrilinos_add_test(fortrilinos_hl/test_eigen_handle.cpp
    ${_needs_anasazi} ${_link})
  fortrilinos_add_test(fortrilinos_hl/test_eigen_handle_preset.cpp
    ${_needs_anasazi} ${_link})
  fortrilinos_add_test(fortrilinos_hl/test_nox.cpp
    ${_link})
  fortrilinos_add_test(fortrilinos_hl/test_solver_handle.cpp
    ${_link})

  fortrilinos_add_test(fortrilinos_hl/test_eigen_handle.F90 ${_needs_anasazi})
  fortrilinos_add_test(fortrilinos_hl/test_eigen_unit.F90 ${_needs_anasazi})
  fortrilinos_add_test(fortrilinos_hl/test_nox.F90)
  fortrilinos_add_test(fortrilinos_hl/test_futility_reproducer.F90)
  fortrilinos_add_test(fortrilinos_hl/test_solver_handle.F90)
  fortrilinos_add_test(fortrilinos_hl/test_solver_unit.F90)
endif()

#---------------------------------------------------------------------------##
