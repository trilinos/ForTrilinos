## Copyright 2017-2018, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE
find_package(Sphinx REQUIRED)

#-----------------------------------------------------------------------------#

# Define a configure file with version and configuration info
if(ForTrilinos_VERSION STREQUAL ForTrilinos_VERSION_STRING)
  set(ForTrilinos_SHORT_VERSION "${ForTrilinos_VERSION}")
elseif("${ForTrilinos_VERSION_STRING}" MATCHES
    "^(.*-rc[.-][0-9]+)(\\.[0-9]+\\+[0-9a-f]+)$")
  # Use release candidate, stripping off stuff since the tag
  set(ForTrilinos_SHORT_VERSION "${CMAKE_MATCH_1}")
else()
  # Development version
  set(_patch "${PROJECT_VERSION_PATCH}")
  if(NOT ForTrilinos_VERSION_STRING MATCHES "-dev")
    # Before the next release
    math(EXPR _patch "${PROJECT_VERSION_PATCH} + 1")
  endif()
  set(ForTrilinos_SHORT_VERSION
    "${PROJECT_VERSION_MAJOR}.${PROJECT_VERSION_MINOR}.${_patch}-dev"
  )
endif()
configure_file("config.json.in" "config.json" @ONLY)

#-----------------------------------------------------------------------------#
add_custom_target(doc)
add_custom_command(TARGET doc
  VERBATIM COMMAND
    "${CMAKE_COMMAND}" -E env
    "CMAKE_CURRENT_BINARY_DIR=${CMAKE_CURRENT_BINARY_DIR}"
  "${SPHINX_EXECUTABLE}" -q
    -d "${CMAKE_CURRENT_BINARY_DIR}/doctrees"
    -b html
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/html"
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  COMMENT "Building HTML documentation with Sphinx"
  DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/conf.py"
  BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
)


install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html/"
  DESTINATION "${CMAKE_INSTALL_DOCDIR}"
)
