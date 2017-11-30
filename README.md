ForTrilinos
===========

[![Build Status](https://cloud.cees.ornl.gov/jenkins-ci/buildStatus/icon?job=ForTrilinos-continuous)](https://cloud.cees.ornl.gov/jenkins-ci/job/ForTrilinos-continuous)
[![Stories in Ready](https://badge.waffle.io/Trilinos/ForTrilinos.svg?label=ready&title=Ready)](http://waffle.io/Trilinos/ForTrilinos)
[![Documentation Status](http://readthedocs.org/projects/fortrilinos/badge/?version=latest)](http://fortrilinos.readthedocs.io/en/latest/?badge=latest)
[![codecov](https://codecov.io/gh/trilinos/ForTrilinos/branch/develop/graph/badge.svg)](https://codecov.io/gh/trilinos/ForTrilinos/branch/develop)

[ForTrilinos](http://trilinos.org/packages/fortrilinos) is an integrated part of the [Trilinos](http://trilinos.org) organization  and provides object-oriented Fortran interfaces to Trilinos C++ packages.

**This code contains the new implementation of the ForTrilinos interfaces using SWIG code generation. At the moment, it is not ready for public consumption as the code and interfaces may change rapidly and backwards compatibility is not guaranteed.**

**The stable version of the original ForTrilinos code developed prior to 2012 is available in [`master`](https://github.com/trilinos/ForTrilinos/tree/master) branch. It is no longer developed or maintained.**

Documentation
-------------

* [Website](http://trilinos.org/packages/fortrilinos)

* [Documentation](http://fortrilinos.readthedocs.org)

Obtaining ForTrilinos
---------------------

ForTrilinos is obtained by cloning the public repository from [GitHub](https://github.com/):

```sh
git clone https://github.com/trilinos/ForTrilinos
```

Installing ForTrilinos
----------------------

Installing ForTrilinos requires:

- Cloning the [Trilinos repository](https://github.com/trilinos/Trilinos) and
  checking out the correct branch:

  ```sh
  git clone git clone https://github.com/trilinos/Trilinos.git
  git checkout trilinos-release-12.12-1
  export TRILINOS_DIR=`pwd`/Trilinos
  ```

- Cloning the [ForTrilinos repository](https://github.com/trilinos/ForTrilinos):

  ```sh
  git clone git clone https://github.com/trilinos/ForTrilinos.git
  export FORTRILINOS_DIR=`pwd`/ForTrilinos
  ```

- Symbolically linking the clone of ForTrilinos to the Trilinos `packages`
  directory:

  ```sh
  ln -s $FORTRILINOS_DIR $TRILINOS_DIR/packages
  ```

- Installing Trilinos with the ForTrilinos package enabled.  See the
  [Installation](http://fortrilinos.readthedocs.io/en/latest/install.html#installation) section of the [ForTrilinos documentation](http://fortrilinos.readthedocs.io/en/latest/index.html) for details on installing Trilinos with ForTrilinos enabled.

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
