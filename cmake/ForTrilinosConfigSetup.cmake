#----------------------------------*-CMake-*----------------------------------#
# Copyright 2020 UT-Battelle, LLC and ForTrilinos developers.
# SPDX-License-Identifier: BSD-3-Clause
# License-Filename: LICENSE
#[=======================================================================[.rst:

ForTrilinosConfigSetup
----------------------

.. command:: fortrilinos_configure_export

  Generate the configure file.

#]=======================================================================]

include(CMakePackageConfigHelpers)

#-----------------------------------------------------------------------------#
# fortrilinos_configure_export
#-----------------------------------------------------------------------------#

macro(fortrilinos_configure_export)
  # Install config files
  install(DIRECTORY "${ForTrilinos_HEADER_CONFIG_DIRECTORY}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
    FILES_MATCHING PATTERN "*.h"
  )
  # Install fortran modules
  install(DIRECTORY "${ForTrilinos_Fortran_MODULE_DIRECTORY}/"
    DESTINATION "${CMAKE_INSTALL_INCLUDEDIR}"
  )

  # Install cmake files
  set(_cmake_files
    "${PROJECT_SOURCE_DIR}/cmake/FindTrilinos.cmake"
  )
  install(FILES ${_cmake_files}
    DESTINATION "${ForTrilinos_INSTALL_CMAKECONFIGDIR}"
  )

  # Export all cache variables that start with ForTrilinos_
  set(ForTrilinos_EXPORT_VARIABLES)
  get_directory_property(_cachevar_keys CACHE_VARIABLES)
  foreach(_key IN LISTS _cachevar_keys)
    if(_key MATCHES "^ForTrilinos_")
      list(APPEND ForTrilinos_EXPORT_VARIABLES "set(${_key} \"${${_key}}\")")
    endif()
  endforeach()

  # Add other cache variables, prefixed with ForTrilinos_
  foreach(_key BUILD_SHARED_LIBS Trilinos_VERSION)
    list(APPEND ForTrilinos_EXPORT_VARIABLES "set(ForTrilinos_${_key} \"${${_key}}\")")
  endforeach()

  # Add hints for TPLs
  foreach(_key MPIEXEC_EXECUTABLE Trilinos_DIR)
    set(_val "${${_key}}")
    if(_val)
      list(APPEND ForTrilinos_EXPORT_VARIABLES "set(${_key} \"${_val}\")")
    endif()
  endforeach()

  list(JOIN ForTrilinos_EXPORT_VARIABLES "\n" ForTrilinos_EXPORT_VARIABLES)

  # Install 'ForTrilinosTargets.cmake', included by ForTrilinosConfig.cmake
  install(EXPORT fortrilinos-targets
    FILE ForTrilinosTargets.cmake
    NAMESPACE ${ForTrilinos_NAMESPACE}
    DESTINATION "${ForTrilinos_INSTALL_CMAKECONFIGDIR}"
  )

  # Generate the file needed by downstream "find_package(ForTrilinos)"
  configure_package_config_file(
    "${PROJECT_SOURCE_DIR}/cmake/ForTrilinosConfig.cmake.in"
    "${ForTrilinos_CMAKE_CONFIG_DIRECTORY}/ForTrilinosConfig.cmake"
    INSTALL_DESTINATION "${ForTrilinos_INSTALL_CMAKECONFIGDIR}"
  )

  # Export version info
  write_basic_package_version_file(
    "${ForTrilinos_CMAKE_CONFIG_DIRECTORY}/ForTrilinosConfigVersion.cmake"
    COMPATIBILITY ExactVersion
  )

  # Install the config and version files
  install(FILES
    "${ForTrilinos_CMAKE_CONFIG_DIRECTORY}/ForTrilinosConfig.cmake"
    "${ForTrilinos_CMAKE_CONFIG_DIRECTORY}/ForTrilinosConfigVersion.cmake"
    DESTINATION ${ForTrilinos_INSTALL_CMAKECONFIGDIR}
  )
endmacro()

#-----------------------------------------------------------------------------#
