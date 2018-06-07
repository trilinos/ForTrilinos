## Copyright 2017-2018, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE
FIND_PACKAGE(Sphinx REQUIRED)

ADD_CUSTOM_TARGET(docs ALL
    ${SPHINX_EXECUTABLE}
        -q -b html
        "${CMAKE_CURRENT_SOURCE_DIR}/source"
        "${CMAKE_CURRENT_BINARY_DIR}/html"
    COMMENT "Building HTML documentation with Sphinx"
)

INSTALL(DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/html"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/docs"
)
