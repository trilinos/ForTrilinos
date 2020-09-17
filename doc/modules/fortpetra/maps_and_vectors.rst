Maps And Vectors
================

Overview
--------

This example demonstrates the following:

- how to create a ``TpetraMap`` object that describes the distribution of other
  Tpetra object entries over processes; and
- how to create the simplest kind of Tpetra linear algebra object: a Vector,
  whose entries are distributed over the process(es) in a communicator

This is a ForTrilinos implementation of the code example associated with Tpetra
`Lesson 02: Map and Vector <https://trilinos.org/docs/dev/packages/tpetra/doc/html/Tpetra_Lesson02.html>`_, which should be consulted for background and theory.

Code example
------------

The following example initializes Tpetra (MPI) then creates two distributed
Tpetra Maps and some Tpetra Vectors, and does a few computations with the
vectors

.. note::

   In the code example below, it is assumed the the following
   ForTrilinos modules are ``use``\ d

   - ``forteuchos``: provides interfaces to the Trilinos' Teuchos package
   - ``fortpetra``: provides interfaces to the Trilinos' Tpetra package
   - ``fortpetra_types``: provides named constants of type default integer for
       Tpetra data types that can be used as the ``kind`` type parameters


.. only:: builder_html

   The source code can be downloaded from
   :download:`here <../../../example/MapsAndVectors.F90>`

.. literalinclude:: ../../../example/MapsAndVectors.F90
   :language: fortran
