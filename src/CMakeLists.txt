
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

SET_AND_INC_DIRS(DIR ${CMAKE_CURRENT_SOURCE_DIR}/epetra)
#APPEND_GLOB(HEADERS ${DIR}/*.hpp ${DIR}/*.h)
APPEND_GLOB(SOURCES ${DIR}/*.f90)

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
  ForTrilinos_enums.f90
  ForTrilinos_utils.f90
  )

#
# C) Define the targets for package's library(s)
#

SET_SOURCE_FILES_PROPERTIES( ${SOURCES} PROPERTIES LANGUAGE Fortran )

PACKAGE_ADD_LIBRARY(
  fortrilinos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

SET_TARGET_PROPERTIES(fortrilinos PROPERTIES LINKER_LANGUAGE Fortran)
