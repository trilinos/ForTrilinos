## Add each new test subdirectory name to the list at the bottom and that's it.

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

## EDIT BLOCK BELOW ONLY
## UNITTEST_* must be at the TOP of the list
ADD_SUBDIRECTORIES(
  UNITTEST_PARSER
  UNITTEST_LAUNCHER
  Epetra_Comm
  Epetra_MpiComm
  Epetra_BlockMap
  Epetra_SerialComm
  )

ENDIF()

######################################################################
######################################################################

