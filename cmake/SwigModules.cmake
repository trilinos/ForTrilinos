set(UseSWIG_MODULE_VERSION 2)
if (CMAKE_VERSION VERSION_LESS 3.20)
  # TODO: This is until Fortran support gets added to the upstream cmake script
  include(UseSWIGFortran)
else()
  cmake_policy(SET CMP0078 "NEW")
  cmake_policy(SET CMP0086 "NEW")
  include(UseSWIG)
endif()

##---------------------------------------------------------------------------##
## ADDING SWIG MODULES
##---------------------------------------------------------------------------##
# MAKE_SWIG_FORTRAN(
#   MODULE module
#   [C]
#   [SOURCE src.i]
#   [DEPLIBS lib1 [lib2 ...]]
#   [DEPMODULES module1 [module2 ...]]
#   [EXTRASRC file1 [file2 ...]]
#   )
#
# Create a SWIG-generated python module and shared object.
#
# Add the [C] flag if the input is to be treated as a C wrapper rather than C++.
#
# The MODULE argument is the name of the resulting module file. By default it
# assumes the name "module.i", but that can be overriden with the SOURCE
# argument.
#
# All libraries in DEPLIBS will be linked against each target.
#
# All other module names in DEPMODULES will be separately added as dependencies
# without being linked. This is used for TriBITS, which understands the
# inter-library linkages but not SWIG modules.
#
# The EXTRASRC argument allows additional sources to be compiled into the SWIG
# module target.

function(MAKE_SWIG_FORTRAN)
  cmake_parse_arguments(PARSE "C" "MODULE;LANGUAGE;SOURCE"
      "DEPLIBS;DEPMODULES;EXTRASRC" ${ARGN})

  if (NOT PARSE_MODULE)
    message(SEND_ERROR "Cannot call MAKE_SWIG without MODULE")
  endif()
  set(PARSE_MODULE ${PARSE_MODULE})

  if (PARSE_SOURCE)
    set(SRC_FILE "${PARSE_SOURCE}")
  else()
    set(SRC_FILE "${PARSE_MODULE}.i")
  endif()

  # Let SWIG know that we're compiling C++ files, and what the module is
  set_source_files_properties(${SRC_FILE} PROPERTIES
    SWIG_MODULE_NAME ${PARSE_MODULE})

  if (NOT PARSE_C)
    set_source_files_properties(${SRC_FILE} PROPERTIES
      CPLUSPLUS TRUE)
  endif()

  # UseSWIG with "NEW" policy creates a target with the same name as the first
  # argument to swig_add_library.
  set(LIBRARY_NAME ${PARSE_MODULE})
  swig_add_library(${LIBRARY_NAME}
    LANGUAGE Fortran
    TYPE USE_BUILD_SHARED_LIBS
    SOURCES ${SRC_FILE} ${PARSE_EXTRASRC})

  # Link against other dependent libraries
  if (PARSE_DEPLIBS)
    target_link_libraries(${LIBRARY_NAME} ${PARSE_DEPLIBS})
  endif()
  foreach(DEPMODULE ${PARSE_DEPMODULES})
    add_dependencies(${LIBRARY_NAME} ${DEPMODULE})
  endforeach()
  set(LINK_LIBS)
# XXX From TribitsLibraryMacros
  TRIBITS_SORT_AND_APPEND_PACKAGE_INCLUDE_AND_LINK_DIRS_AND_LIBS(
    ${PACKAGE_NAME}  LIB  LINK_LIBS)

  TRIBITS_SORT_AND_APPEND_TPL_INCLUDE_AND_LINK_DIRS_AND_LIBS(
    ${PACKAGE_NAME}  LIB  LINK_LIBS)
# XXX END

  # Set properties on the target.
  get_target_property(INCL_DIR ${LIBRARY_NAME} INCLUDE_DIRECTORIES)
  if (NOT INCL_DIR)
    set(INCL_DIR)
  endif()
  list(APPEND INCL_DIR ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR})
  list(REMOVE_DUPLICATES INCL_DIR)
  set_target_properties(${LIBRARY_NAME} PROPERTIES
    SWIG_INCLUDE_DIRECTORIES "${INCL_DIR}"
    COMPILE_FLAGS "${SWIG_CXX_FLAGS}")
  target_link_libraries(${LIBRARY_NAME} ${LINK_LIBS})

# XXX From TribitsLibraryMacros
  # Add to tribits library list
  SET_PROPERTY(
    TARGET ${LIBRARY_NAME}
    APPEND PROPERTY
    LABELS ${PACKAGE_NAME}Libs ${PARENT_PACKAGE_NAME}Libs
    )
  PREPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_LIB_TARGETS ${LIBRARY_NAME})
  PREPEND_GLOBAL_SET(${PARENT_PACKAGE_NAME}_ALL_TARGETS ${LIBRARY_NAME})
  INSTALL(
    TARGETS ${LIBRARY_NAME}
    EXPORT ${PACKAGE_NAME}
    RUNTIME DESTINATION "${${PROJECT_NAME}_INSTALL_RUNTIME_DIR}"
    LIBRARY DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
    ARCHIVE DESTINATION "${${PROJECT_NAME}_INSTALL_LIB_DIR}"
    COMPONENT ${PACKAGE_NAME}
    )
  PREPEND_GLOBAL_SET(${PACKAGE_NAME}_INCLUDE_DIRS  ${INCL_DIR})
  PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARY_DIRS  ${INCL_DIR})
  PREPEND_GLOBAL_SET(${PACKAGE_NAME}_LIBRARIES  ${LIBRARY_NAME})

  GLOBAL_SET(${PACKAGE_NAME}_HAS_NATIVE_LIBRARIES TRUE)
# XXX END TriBits

  # Install the built target module
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${PARSE_MODULE}.mod
    DESTINATION include)
endfunction()


