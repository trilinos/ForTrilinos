
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL")
    SET(CMAKE_Fortran_FLAGS "${CMake_Fortran_FLAGS} -qxlf2003=polymorphic -qxlf2003=autorealloc -WF,-C! -tchbI -B/opt/local/hotfix/091023/exe/")
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

SET_AND_INC_DIRS(DIREPETRA ${CMAKE_CURRENT_SOURCE_DIR}/epetra)
SET_AND_INC_DIRS(DIRTEUCHOS ${CMAKE_CURRENT_SOURCE_DIR}/teuchos)
SET_AND_INC_DIRS(DIRAMESOS ${CMAKE_CURRENT_SOURCE_DIR}/amesos)
SET_AND_INC_DIRS(DIRAZTECOO ${CMAKE_CURRENT_SOURCE_DIR}/aztecoo)
SET_AND_INC_DIRS(DIRGALERI ${CMAKE_CURRENT_SOURCE_DIR}/galeri)
SET_AND_INC_DIRS(DIRIFPACK ${CMAKE_CURRENT_SOURCE_DIR}/ifpack)

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

# Set sources for link-layer build

APPEND_SET(SOURCES
  ForTrilinos_enum_wrappers.F90
  ForTrilinos_enums.F90
  ForTrilinos_utils.F90
  ${DIREPETRA}/forepetra.F90
  ${DIRTEUCHOS}/forteuchos.F90
  )

IF (${PACKAGE_NAME}_ENABLE_Amesos)

  APPEND_SET(SOURCES
    ${DIRAMESOS}/foramesos.F90
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_AztecOO)

  APPEND_SET(SOURCES
    ${DIRAZTECOO}/foraztecoo.F90
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Galeri)

  APPEND_SET(SOURCES
    ${DIRGALERI}/forgaleri.F90
    )

ENDIF()

IF (${PACKAGE_NAME}_ENABLE_Ifpack)

  APPEND_SET(SOURCES
    ${DIRIFPACK}/forifpack.F90
    )

ENDIF()

# Set additional sources for object-oriented build

IF(${PACKAGE_NAME}_ENABLE_EXTENDED)

  APPEND_SET(SOURCES
    ForTrilinos_hermetic.F90
    ForTrilinos_universal.F90
    )

  APPEND_GLOB(SOURCES ${DIREPETRA}/FEpetra*.F90
                      ${DIRTEUCHOS}/FTeuchos*.F90)

  IF (${PACKAGE_NAME}_ENABLE_Amesos)

    APPEND_GLOB(SOURCES ${DIRAMESOS}/FAmesos*.F90)

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_AztecOO)

    APPEND_GLOB(SOURCES ${DIRAZTECOO}/FAztecOO*.F90)

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Galeri)

    APPEND_GLOB(SOURCES ${DIRAZTECOO}/FGaleri*.F90)

  ENDIF()

  IF (${PACKAGE_NAME}_ENABLE_Ifpack)

    APPEND_GLOB(SOURCES ${DIRAZTECOO}/FIfpack*.F90)

  ENDIF()

ENDIF()

#
# C) Define the targets for package's library(s)
#

PACKAGE_ADD_LIBRARY(
  fortrilinos
  HEADERS ${HEADERS}
  SOURCES ${SOURCES}
  )

SET_TARGET_PROPERTIES(fortrilinos PROPERTIES LINKER_LANGUAGE Fortran)
