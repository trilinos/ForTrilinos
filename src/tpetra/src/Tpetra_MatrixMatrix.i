/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%import <std_string.i>

%{
#include "Teuchos_RCP.hpp"
#include "TpetraExt_MatrixMatrix.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================
%ignore Tpetra::MatrixMatrix::add;          // Same interface as 2-arg Add leads to conflict

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::MatrixMatrix::Jacobi;       // needs Tpetra::Vector

%include "TpetraExt_MatrixMatrix_decl.hpp"

%template(TpetraMatrixMatrixMultiply) Tpetra::MatrixMatrix::Multiply<SC,LO,GO,NO>;
%template(TpetraMatrixMatrixAdd)      Tpetra::MatrixMatrix::Add<SC,LO,GO,NO>;
