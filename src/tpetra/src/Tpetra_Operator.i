/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"
%import <std_string.i>
%import <Teuchos_Comm.i>

%{
#include "Teuchos_RCP.hpp"
#include "Tpetra_Operator.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================

// =======================================================================
// Postpone temporarily
// =======================================================================

%include "Tpetra_Operator.hpp"

%teuchos_rcp(Tpetra::Operator<SC,LO,GO,NO>)
%template(TpetraOperator) Tpetra::Operator<SC,LO,GO,NO>;
