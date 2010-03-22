
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  verySimple
  verySimpleObjectOriented
  FEpetra_power_method
  )

IF(${PACKAGE_NAME}_ENABLE_EXTENDED)
  ADD_SUBDIRECTORIES(
    )
ENDIF()
