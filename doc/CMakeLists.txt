## Copyright 2017-2018, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE

configure_file("version.json.in" "version.json" @ONLY)
add_custom_target(doc)
add_custom_command(TARGET doc
  COMMAND "${CMAKE_COMMAND}" -E env "BUILD_DIR=${CMAKE_CURRENT_BINARY_DIR}"
  "${SPHINX_EXECUTABLE}" -q
    -d "${CMAKE_CURRENT_BINARY_DIR}/doctrees"
    -b html
    "${CMAKE_CURRENT_SOURCE_DIR}"
    "${CMAKE_CURRENT_BINARY_DIR}/html"
  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
  COMMENT "Building HTML documentation with Sphinx"
  BYPRODUCTS "${CMAKE_CURRENT_BINARY_DIR}/html/index.html"
)

install(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html/"
  DESTINATION "${CMAKE_INSTALL_DOCDIR}"
)
