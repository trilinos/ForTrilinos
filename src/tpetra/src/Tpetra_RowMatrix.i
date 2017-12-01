/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"
%import <Teuchos_Comm.i>

%{
#include "Tpetra_RowMatrix.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================
%ignore Tpetra::RowMatrix::getNode;
// Ignore implementation details
%ignore Tpetra::RowMatrix::pack;

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::RowMatrix::getFrobeniusNorm;     // needs to understand that Tpetra::MultiVector<SC,LO,GO,NO,false> and Tpetra::MultiVector<SC,LO,GO,NO> are the same thing
%ignore Tpetra::RowMatrix::getGraph;             // needs Tpetra::RowGraph
%ignore Tpetra::RowMatrix::getLocalDiagCopy;     // needs Tpetra::Vector
%ignore Tpetra::RowMatrix::leftScale;            // needs Tpetra::Vector
%ignore Tpetra::RowMatrix::rightScale;           // needs Tpetra::Vector


%teuchos_rcp(Tpetra::RowMatrix<SC,LO,GO,NO>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_RowMatrix_decl.hpp"

%template(TpetraRowMatrix) Tpetra::RowMatrix<SC,LO,GO,NO>;
