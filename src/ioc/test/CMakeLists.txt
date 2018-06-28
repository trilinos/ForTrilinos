## Copyright 2017-2018, UT-Battelle, LLC
##
## SPDX-License-Identifier: BSD-3-Clause
## License-Filename: LICENSE
TRIBITS_COPY_FILES_TO_BINARY_DIR(IOC_cp
  SOURCE_FILES params.xml
 )

TRIBITS_ADD_EXECUTABLE(
  test_tpetra_fem_ioc_cpp
  SOURCES test_tpetra_fem_ioc.cpp
  COMM serial mpi
)

