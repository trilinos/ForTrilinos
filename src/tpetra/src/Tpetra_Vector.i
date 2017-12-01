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
#include "Tpetra_Vector.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::Vector::assign;
%ignore Tpetra::Vector::Vector(const Teuchos::RCP<const map_type> &map, const dual_view_type &view);    // needs Kokkos::DualView
%ignore Tpetra::Vector::Vector(const Teuchos::RCP<const map_type> &map, const dual_view_type &view, const dual_view_type &origView);    // needs Kokkos::DualView
%ignore Tpetra::Vector::getData;                // needs Teuchos::ArrayRCP
%ignore Tpetra::Vector::getDataNonConst;        // needs Teuchos::ArrayRCP
%ignore Tpetra::Vector::describe;               // needs Teuchos::FancyOStream

// =======================================================================
// Fix ±1 issues
// =======================================================================
%typemap(in)  int myRow %{$1 = *$input - 1;%}

%teuchos_rcp(Tpetra::Vector<SC,LO,GO,NO,false>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_Vector_decl.hpp"

%template(TpetraVector) Tpetra::Vector<SC,LO,GO,NO,false>;
