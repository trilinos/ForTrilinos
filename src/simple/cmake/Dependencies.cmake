TRIBITS_PACKAGE_DEFINE_DEPENDENCIES(
  LIB_REQUIRED_PACKAGES
    Stratimikos
    Teuchos
    Thyra
    ThyraTpetraAdapters
    Tpetra
    ForTrilinosTeuchos
  LIB_OPTIONAL_PACKAGES
    Ifpack2
    MueLu
  TEST_REQUIRED_PACKAGES
    Belos # Needed to pass, not to build
  TEST_OPTIONAL_PACKAGES
    Belos
    Ifpack2
    MueLu
)
