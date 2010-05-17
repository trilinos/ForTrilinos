## Add each new test subdirectory name to the list at the bottom and that's it.

INCLUDE(PackageAddTest)
INCLUDE(AddSubdirectories)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

######################################################################

IF(FORTRILINOS_COMPILER_SUPPORTS_LONG_LINES)
  SET(FORTRILINOS_LINE_BREAKER
    -DFT_SUPPORTS_LONG_LINES
    )
ELSE()
  MESSAGE("*** ForTrilinos unit tests will be skipped! ***")
  SET(FORTRILINOS_SKIP_UNITTESTS 1)
#    SET(FORTRILINOS_LINE_BREAKER
#      | sed 's/@@/\\n/g'
#      )
ENDIF()

######################################################################
######################################################################

IF(NOT FORTRILINOS_SKIP_UNITTESTS)

######################################################################

PACKAGE_ADD_EXECUTABLE(unittest_launcher SOURCES unittest_launcher.cpp)

SET_TARGET_PROPERTIES(
  ForTrilinos_unittest_launcher PROPERTIES
  LINKER_LANGUAGE CXX
  )

######################################################################

SET(ALL_MACRO_FILES
  ${CMAKE_CURRENT_SOURCE_DIR}/all_build_macros.h
  ${CMAKE_CURRENT_SOURCE_DIR}/cxx_driver_macros.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_driver_macros.h
  ${CMAKE_CURRENT_SOURCE_DIR}/fortran_impl_macros.h
  ${CMAKE_CURRENT_SOURCE_DIR}/runtime_macros.h
  ${CMAKE_CURRENT_SOURCE_DIR}/breaklines.h
  )

######################################################################

## EDIT BLOCK BELOW ONLY
ADD_SUBDIRECTORIES(
  Epetra_Comm
  )

######################################################################

ENDIF()

######################################################################
######################################################################

