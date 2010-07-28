
INCLUDE(PackageMacros)

PACKAGE(ForTrilinos Fortran)

IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qfree=f90 -qxlf2003=polymorphic -qxlf2003=autorealloc -WF,-C!")
#    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -qxlf2003=polymorphic -qxlf2003=autorealloc -qxlf90=autodealloc -WF,-C!")
    SET(ForTrilinos_UT_CPP_FLAGS -E -P -qnoppline)
    SET(ForTrilinos_UT_CPP ${CMAKE_C_COMPILER})
ELSE()
    SET(ForTrilinos_UT_CPP_FLAGS "")
    SET(ForTrilinos_UT_CPP cpp)
ENDIF()

# Unit testing requires support for long lines
# Enabled by default for: "NAG", "XL", Intel
# Limited support: "PGI" (limited to ~250)
# Flag needed for: "GNU"
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "GNU")
    SET(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -ffree-line-length-none")
ENDIF()
# If not one of these compilers, will need to break lines using sed
IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL|NAG|Intel|GNU")
    SET(FORTRILINOS_COMPILER_SUPPORTS_LONG_LINES 1)
ENDIF()

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_POSTPROCESS()
