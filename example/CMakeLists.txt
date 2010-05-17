
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  verySimple
  )

IF(${PACKAGE_NAME}_ENABLE_EXTENDED)
  ADD_SUBDIRECTORIES(
  verySimpleObjectOriented
  FEpetra_power_method
    )
ENDIF()
