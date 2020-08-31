/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%{
#include <Tpetra_Vector.hpp>
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

%include "Tpetra_Vector_decl.hpp"

%teuchos_rcp(Tpetra::Vector<SC,LO,GO,NO>)
%template(TpetraVector) Tpetra::Vector<SC,LO,GO,NO>;
