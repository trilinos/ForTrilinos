ForTrilinos
===========

[![Build Status](https://cloud.cees.ornl.gov/jenkins-ci/buildStatus/icon?job=ForTrilinos-continuous)](https://cloud.cees.ornl.gov/jenkins-ci/job/ForTrilinos-continuous)
[![Documentation Status](http://readthedocs.org/projects/fortrilinos/badge/?version=latest)](http://fortrilinos.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/trilinos/ForTrilinos/branch/develop/graph/badge.svg)](https://codecov.io/gh/trilinos/ForTrilinos/branch/develop)

[ForTrilinos](http://trilinos.org/packages/fortrilinos) is a part of the [Trilinos](http://trilinos.org) project and provides object-oriented Fortran interfaces to Trilinos C++ packages.

This is the new effort to provide Fortran interfaces to Trilinos through
automatic code generation using SWIG. The previous effort (ca. 2008-2012) can
be obtained by downloading Trilinos releases prior to 12.12.

Provided functionality
----------------------
ForTrilinos provides Fortran interfaces for the following capabilities:
- Parameter lists and XML parsers (through Teuchos);
- Distributed linear algebra object including sparse graphs, sparse matrices, and dense vectors (through Tpetra);
- Linear solvers and preconditioners (through Stratimikos, Ifpack2, Belos, MueLu);
- Eigen solvers (through Anasazi).

Documentation
-------------

* [Website](http://trilinos.org/packages/fortrilinos)

* [Documentation](http://fortrilinos.readthedocs.org)

Installing ForTrilinos
----------------------

Please consult the documentation available [here](https://fortrilinos.readthedocs.io/en/latest/install.html).

Questions, Bug Reporting, and Issue Tracking
--------------------------------------------

Questions, bug reporting and issue tracking are provided by GitHub. Please
report all bugs by creating a new issue with the bug tag. You can ask
questions by creating a new issue with the question tag.

Contributing
------------
We encourage you to contribute to ForTrilinos! Please check out the
[guidelines](CONTRIBUTING.md) about how to proceed.

License
-------
ForTrilinos is licensed under a BSD license.
