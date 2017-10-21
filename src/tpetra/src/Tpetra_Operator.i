// Dependencies
%include "Teuchos_RCP.i"
%import <std_string.i>
%import <Teuchos_ArrayView.i>
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

%teuchos_rcp(Tpetra::Operator<SC,LO,GO,NO>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_Operator.hpp"

%template(TpetraOperator) Tpetra::Operator<SC,LO,GO,NO>;
