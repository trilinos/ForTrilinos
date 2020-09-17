Modules
=======

ForTrilinos includes several Fortran modules (``forteuchos``,
``forbelos``, ``fortpetra``) that are thin layers built on top of
Trilinos packages. The ``fortrilinos_hl`` package is a high-level set of
wrappers that exposes linear, nonlinear, and eigenvalue solvers.

Different Trilinos packages are required for different levels of
functionality:

=========== =============================
Package     Module
=========== =============================
Teuchos     ``forteuchos``
Belos       ``forbelos``
Tpetra      ``fortpetra``
Anasazi     ``fortrilinos_hl``
NOX         ``fortrilinos_hl``
Stratimikos ``fortrilinos_hl``
Thyra       ``fortrilinos_hl``
Amesos2     ``fortrilinos_hl`` (optional)
Ifpack2     ``fortrilinos_hl`` (optional)
MueLu       ``fortrilinos_hl`` (optional)
=========== =============================

.. toctree::
   :maxdepth: 1

   modules/fortpetra.rst
