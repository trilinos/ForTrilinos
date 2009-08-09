
INCLUDE(PackageMacros)

PACKAGE(ForTrilinos Fortran)

#IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL")
#  SET(CMAKE_Fortran_COMPILER xlf2003)
#ENDIF()

IF("${CMAKE_Fortran_COMPILER_ID}" MATCHES "VisualAge|XL")
    SET(CMAKE_Fortran_FLAGS "${CMake_Fortran_FLAGS} -qxlf2003=polymorphic -qxlf2003=autorealloc -WF,-C!")
ENDIF()

ADD_SUBDIRECTORY(src)

PACKAGE_ADD_TEST_DIRECTORIES(test)

#PACKAGE_ADD_EXAMPLE_DIRECTORIES(example)

PACKAGE_POSTPROCESS()
