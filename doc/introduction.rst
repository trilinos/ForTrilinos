************
Introduction
************

ForTrilinos is an open-source software library providing object-oriented
Fortran interfaces to the Trilinos parallel numerical software collection. Its
objective is to facilitate highly parallel Fortran application development by
exposing an idiomatic library numerical algorithms.

From about 2008 to 2012, Trilinos included a set of hand-rolled wrappers to the
first-generation solver packages Epetra, Amesos, Ifpack, and AztecOO
:cite:`OldForTrilinos`. Support for these older packages and the Fortran
interfaces faded away, and at the inception of the Exascale Computing Project
:cite:`Ecp` several Fortran codes intended to utilize the new heterogeneous
architecture capabilities in Trilinos, motivating a new version targeting
the newer solver libraries and focusing on maintainability.

This need for reducing the maintenance burden motivated development of a
Fortran extension to SWIG, a tool for generating inter-language bindings to C++
libraries :cite:`SwigFortran`.  The new SWIG-based implementation of
ForTrilinos, numbered starting at version 2 to acknowledge the prior set
of Fortran bindings, started in late 2016 and is being actively maintained as
part of ECP. As of version 2.0.0, the documented ForTrilinos module interfaces
are stable and available for public release.
