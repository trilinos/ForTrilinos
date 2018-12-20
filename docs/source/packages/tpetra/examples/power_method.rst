Power Method for Finding the Largest Eigenvalue
===============================================

Overview
--------

This example demonstrates the following:

- how to construct a ``TpetraCrsMatrix`` (a distributed sparse matrix);
- how to modify the entries of a previously constructed ``TpetraCrsMatrix``; and
- how to use ``TpetraCrsMatrix`` and ``TpetraMultiVector`` to implement a simple
  iterative eigensolver (the power method)

This example is a ForTrilinos implementation of the code example associated with
Tpetra `Lesson 03: Power Method <https://trilinos.org/docs/dev/packages/tpetra/doc/html/Tpetra_Lesson03.html>`_, which should be consulted for background and theory.

Code example
------------

The following code example shows how to fill and compute with a Tpetra sparse
matrix, using the procedure discussed in the text above.

.. note::

   In the code example below, it is assumed the the following
   ForTrilinos modules are ``use``\ d

   - ``forteuchos``: provides interfaces to the Trilinos' Teuchos package
   - ``fortpetra``: provides interfaces to the Trilinos' Tpetra package
   - ``fortpetra_types``: provides named constants of type default integer for
       Tpetra data types that can be used as the ``kind`` type parameters


.. only:: builder_html

   The source code can be downloaded from
   :download:`here <../../../../../src/tpetra/tutorials/PowerMethod.F90>`

.. literalinclude:: ../../../../../src/tpetra/tutorials/PowerMethod.F90
   :language: fortran
