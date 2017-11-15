INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR})

IF (ForTrilinos_ENABLE_DeveloperMode)
  IF (ForTrilinos_EXPORT_SWIG)
    # Generate wrapper and proxy files since MPI is enabled
    SET(CMAKE_SWIG_OUTDIR "${CMAKE_CURRENT_SOURCE_DIR}")
  ENDIF()

  MAKE_SWIG(MODULE forbelos
    SOURCE forbelos.i
    EXTRASRC ${EXTRA_SOURCE}
    )

ELSE()
  # Use pre-generated files
  SET(GEN_SOURCE
    forbelos.F90
    forbelosFORTRAN_wrap.cxx)

  TRIBITS_ADD_LIBRARY(
    forbelos
    SOURCES ${GEN_SOURCE} ${EXTRA_SOURCE})
ENDIF()

INSTALL(FILES ${CMAKE_CURRENT_BINARY_DIR}/forbelos.mod DESTINATION include)
