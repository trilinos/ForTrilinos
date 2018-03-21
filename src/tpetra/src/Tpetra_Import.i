/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"

%{
#include "Tpetra_Import.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================
// Ignore Details namespace

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::Import::Import(const Teuchos::RCP< const map_type > &source, \
        const Teuchos::RCP< const map_type > &target, \
        const Teuchos::RCP< Teuchos::FancyOStream > &out);      // needs Teuchos::FancyOStream
%ignore Tpetra::Import::Import(const Teuchos::RCP< const map_type > &source, \
        const Teuchos::RCP< const map_type > &target, \
        const Teuchos::RCP< Teuchos::FancyOStream > &out, \
        const Teuchos::RCP< Teuchos::ParameterList > &plist);   // needs Teuchos::FancyOStream
%ignore Tpetra::Import::Import (const Teuchos::RCP<const map_type>& source,
        const Teuchos::RCP<const map_type>& target,
        Teuchos::Array<int> & remotePIDs);              // ±1 issue
%ignore Tpetra::Import::Import( \
        const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &source, \
        const Teuchos::RCP< const Map< LocalOrdinal, GlobalOrdinal, Node > > &target, \
        Teuchos::Array< int > &userRemotePIDs, Teuchos::Array< GlobalOrdinal > &remoteGIDs, \
        const Teuchos::ArrayView< const LocalOrdinal > &userExportLIDs, \
        const Teuchos::ArrayView< const int > &userExportPIDs, \
        const bool useRemotePIDs, \
        const Teuchos::RCP< Teuchos::ParameterList > &plist=Teuchos::null, \
        const Teuchos::RCP< Teuchos::FancyOStream > &out=Teuchos::null); // ±1 issue, needs Teuchos::FancyOStream
%ignore Tpetra::Import::getPermuteFromLIDs;     // ±1 issue
%ignore Tpetra::Import::getPermuteToLIDs;       // ±1 issue
%ignore Tpetra::Import::getRemoteLIDs;          // ±1 issue
%ignore Tpetra::Import::getExportLIDs;          // ±1 issue
%ignore Tpetra::Import::getExportPIDs;          // ±1 issue
%ignore Tpetra::Import::getDistributor;         // needs Tpetra::Distributor
%ignore Tpetra::Import::operator=;              // needs operator=
%ignore Tpetra::Import::describe;               // needs Teuchos::FancyOStream
%ignore Tpetra::Import::print;                  // needs std::ostream


%teuchos_rcp(Tpetra::Import<LO,GO,NO>)

%include "Tpetra_Import_decl.hpp"

%template(TpetraImport) Tpetra::Import<LO,GO,NO>;
