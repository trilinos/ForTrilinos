Overview
========

Tpetra_

* is a package within the Trilinos_ ecosystem;
* implements linear algebra objects, including sparse graphs, sparse matrices,
  and dense vectors;
* provides parallel data redistribution facilities for many other Trilinos'
  packages; and
* is “hybrid parallel,” meaning that it uses at least two levels of parallelism.

Within Trilinos, Tpetra can be used mostly as a stand-alone package, with
explicit dependencies on Teuchos_ and Kokkos_. There are adapters allowing the use
of Tpetra operators and multivectors in both the Belos_ linear solver package and
the Anasazi_ eigensolver package.

Hybrid Parallel
---------------

Hybrid parallel refers to the fact that Tpetra uses at least two levels of parallelism:

* MPI_ (the Message Passing Interface) for distributed-memory parallelism; and
* any of various shared-memory parallel programming models within an MPI process.

Tpetra uses the Kokkos_ package’s shared-memory parallel programming model for
data structures and computational kernels.  Kokkos makes it easy to port Tpetra
to new computer architectures and to extend its use of parallel computational
kernels and thread-scalable data structures.  Kokkos supports several different
parallel programming models, including

- OpenMP_
- `POSIX Threads`_ (Pthreads)
- Nvidia_’s CUDA_ programming model for graphics processing units (GPUs)

Features
--------

Tpetra has the following unique features:

* Native support for representing and solving very large graphs, matrices, and
  vectors.  “Very large” means over two billion unknowns or other entities.
* Matrices and vectors may contain many different kinds of data, such as
  floating-point types of different precision, complex-valued types, automatic
  differentiation objects from the Sacado package, or stochastic PDE
  discretization types from the Stokhos package.
* Support for many different shared-memory parallel programming models

Tpetra Support in ForTrilinos
-----------------------------

Tpetra is a large library and its full functionality has not been exposed in
ForTrilinos, though coverage is actively being expanded.  The Tpetra
:ref:`tpetra_examples` provide a good overview of the Tpetra functionality currently
supported.

.. include:: /links.txt
