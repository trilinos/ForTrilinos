
INCLUDE(PackageAddExecutableAndTest)
INCLUDE(AddSubdirectories)

ADD_SUBDIRECTORIES(
  verySimple
  )

IF(${PACKAGE_NAME}_ENABLE_EXTENDED)

  IF (${PACKAGE_NAME}_ENABLE_AztecOO)
    ADD_SUBDIRECTORIES(
      aztecoo
       )
  ENDIF()

  ADD_SUBDIRECTORIES(
   verySimpleObjectOriented
   FEpetra_power_method
    )
ENDIF()
