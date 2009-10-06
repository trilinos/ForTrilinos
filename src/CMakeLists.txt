
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL")
    SET(CMAKE_Fortran_FLAGS "${CMake_Fortran_FLAGS} -qxlf2003=polymorphic -qxlf2003=autorealloc -WF,-C!")
ENDIF()

INCLUDE(PackageLibraryMacros)

#
# A) Package-specific configuration options
#

PACKAGE_CONFIGURE_FILE(${PACKAGE_NAME}_config.h)

#
# B) Define the header and source files (and include directories)
#

SET(HEADERS "")
SET(SOURCES "")

# Set the base Trilinos directory in order to find Trilinos_version.h

ASSERT_DEFINED(PROJECT_SOURCE_DIR)

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR})

APPEND_SET(HEADERS
  ${PROJECT_SOURCE_DIR}/Trilinos_version.h
  )

# Set headers (for IDEs and installing)

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

APPEND_SET(HEADERS
  ${CMAKE_CURRENT_BINARY_DIR}/${PACKAGE_NAME}_config.h
  )

# NOTE: Above, we add the XXX_config.h file so tht MS Visual Studio
# will know to include this file in the dependency analysis

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})

APPEND_SET(HEADERS
  )

# Set sources (add any source files that in the current directory [src])

APPEND_SET(SOURCES
  ForTrilinos_enums.F90
  epetra/forepetra.F90
  teuchos/forteuchos.F90
  amesos/foramesos.F90
  ForTrilinos_utils.F90
#  ForTrilinos_hermetic.F90
#  ForTrilinos_universal.F90
  )

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  fortrilinos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

SET_TARGET_PROPERTIES(fortrilinos PROPERTIES LINKER_LANGUAGE Fortran)
