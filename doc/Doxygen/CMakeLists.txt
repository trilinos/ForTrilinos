# -*- cmake -*-

# @HEADER
# ***********************************************************************
#
#          ForTrilinos: Fortran Interfaces to Trilinos Packages
# ***********************************************************************
# @HEADER

SET(ADDITIONAL_CLEAN_FILES "")

FOREACH(Package ${ForTrilinos_PACKAGES})

  IF(NOT EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${Package}_dox.i)

    IF(ForTrilinos_DOCSTRINGS)

      IF(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_${Package}.in)

        CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile_${Package}.in
          ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile_${Package}
          @ONLY
          )
        EXECUTE_PROCESS(COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/docstrings.py
          --doxygen=${DOXYGEN_EXECUTABLE} ${Package}
          WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

      ENDIF()

    ELSE(ForTrilinos_DOCSTRINGS)

      EXECUTE_PROCESS(COMMAND ${CMAKE_COMMAND} -E touch
        ${CMAKE_CURRENT_BINARY_DIR}/${Package}_dox.i
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})
      MESSAGE(STATUS
        "ForTrilinos.${Package} will not have doxygen-generated docstrings")

    ENDIF(ForTrilinos_DOCSTRINGS)

    IF(EXISTS ${CMAKE_CURRENT_BINARY_DIR}/${Package}_dox.i)
      LIST(APPEND ADDITIONAL_CLEAN_FILES ${Package}_dox.i)
    ENDIF()

  ENDIF()

ENDFOREACH(Package)

# Set this directory's clean files
GET_DIRECTORY_PROPERTY(clean_files ADDITIONAL_MAKE_CLEAN_FILES)
LIST(APPEND            clean_files ${ADDITIONAL_CLEAN_FILES})
LIST(REMOVE_DUPLICATES clean_files)
LIST(REMOVE_ITEM       clean_files "")
SET_DIRECTORY_PROPERTIES(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${clean_files}")
