/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"

%{
#include "Tpetra_Export.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================
// Ignore Details namespace

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::Export::Export(const Teuchos::RCP< const map_type > &source, \
        const Teuchos::RCP< const map_type > &target, \
        const Teuchos::RCP< Teuchos::FancyOStream > &out);      // needs Teuchos::FancyOStream
%ignore Tpetra::Export::Export(const Teuchos::RCP< const map_type > &source, \
        const Teuchos::RCP< const map_type > &target, \
        const Teuchos::RCP< Teuchos::FancyOStream > &out, \
        const Teuchos::RCP< Teuchos::ParameterList > &plist);   // needs Teuchos::FancyOStream
%ignore Tpetra::Export::getPermuteFromLIDs;     // ±1 issue
%ignore Tpetra::Export::getPermuteToLIDs;       // ±1 issue
%ignore Tpetra::Export::getRemoteLIDs;          // ±1 issue
%ignore Tpetra::Export::getExportLIDs;          // ±1 issue
%ignore Tpetra::Export::getExportPIDs;          // ±1 issue
%ignore Tpetra::Export::getDistributor;         // needs Tpetra::Distributor
%ignore Tpetra::Export::operator=;              // needs operator=
%ignore Tpetra::Export::describe;               // needs Teuchos::FancyOStream
%ignore Tpetra::Export::print;                  // needs std::ostream


%teuchos_rcp(Tpetra::Export<LO,GO,NO>)

%include "Tpetra_Export_decl.hpp"

%template(TpetraExport) Tpetra::Export<LO,GO,NO>;
