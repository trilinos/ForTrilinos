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
  TEST_REQUIRED_PACKAGES # Needed to pass, not to build
    Belos
    Ifpack2
  TEST_OPTIONAL_PACKAGES
    MueLu
)
