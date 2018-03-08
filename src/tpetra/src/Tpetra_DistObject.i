/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"

%{
#include "Tpetra_DistObject.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================
%ignore Tpetra::DistObject::DistObject;

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::DistObject::print;                  // needs std::ostream
%ignore Tpetra::DistObject::describe;               // needs Teuchos::FancyOStream

%teuchos_rcp(Tpetra::SrcDistObject);
%include "Tpetra_SrcDistObject.hpp"

%include "Tpetra_DistObject_decl.hpp"
%teuchos_rcp(Tpetra::DistObject<Packet,LO,GO,NO,false>);
%template(TpetraDistObject) Tpetra::DistObject<Packet,LO,GO,NO,false>;
