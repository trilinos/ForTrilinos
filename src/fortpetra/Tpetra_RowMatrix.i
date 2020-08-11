/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%{
#include <Tpetra_RowMatrix.hpp>
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


%include "Tpetra_RowMatrix_decl.hpp"

%teuchos_rcp(Tpetra::RowMatrix<SC,LO,GO,NO>)
%template(TpetraRowMatrix) Tpetra::RowMatrix<SC,LO,GO,NO>;
