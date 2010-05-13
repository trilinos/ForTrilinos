## Add each new test subdirectory name to the list at the bottom and that's it.

INCLUDE(PackageAddTest)
INCLUDE(AddSubdirectories)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  )

######################################################################

PACKAGE_ADD_EXECUTABLE(unittest_launcher SOURCES unittest_launcher.cpp)

SET_TARGET_PROPERTIES(
  ForTrilinos_unittest_launcher PROPERTIES
  LINKER_LANGUAGE CXX
  )

######################################################################

## EDIT BLOCK BELOW ONLY
ADD_SUBDIRECTORIES(
  Epetra_Comm
  )

######################################################################

