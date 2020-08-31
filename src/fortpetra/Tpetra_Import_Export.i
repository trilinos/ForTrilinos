/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%{
#include <Tpetra_Import.hpp>
#include <Tpetra_Export.hpp>
%}

// =======================================================================
// Ignore permanently
// =======================================================================
// Ignore Details namespace

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::Export::Export(const Teuchos::RCP< const map_type > &source,
        const Teuchos::RCP< const map_type > &target,
        const Teuchos::RCP< Teuchos::FancyOStream > &out);      // needs Teuchos::FancyOStream
%ignore Tpetra::Export::Export(const Teuchos::RCP< const map_type > &source,
        const Teuchos::RCP< const map_type > &target,
        const Teuchos::RCP< Teuchos::FancyOStream > &out,
        const Teuchos::RCP< Teuchos::ParameterList > &plist);   // needs Teuchos::FancyOStream
%ignore Tpetra::Export::getPermuteFromLIDs;     // ±1 issue
%ignore Tpetra::Export::getPermuteToLIDs;       // ±1 issue
%ignore Tpetra::Export::getRemoteLIDs;          // ±1 issue
%ignore Tpetra::Export::getExportLIDs;          // ±1 issue
%ignore Tpetra::Export::getExportPIDs;          // ±1 issue
%ignore Tpetra::Export::getDistributor;         // needs Tpetra::Distributor

%ignore Tpetra::Import::Import(const Teuchos::RCP< const map_type > &source,
        const Teuchos::RCP< const map_type > &target,
        const Teuchos::RCP< Teuchos::FancyOStream > &out);      // needs Teuchos::FancyOStream
%ignore Tpetra::Import::Import(const Teuchos::RCP< const map_type > &source,
        const Teuchos::RCP< const map_type > &target,
        const Teuchos::RCP< Teuchos::FancyOStream > &out,
        const Teuchos::RCP< Teuchos::ParameterList > &plist);   // needs Teuchos::FancyOStream
%ignore Tpetra::Import::Import (
        const Teuchos::RCP<const map_type>& source,
        const Teuchos::RCP<const map_type>& target,
        Teuchos::Array<int> & remotePIDs,
        const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::rcp(new Teuchos::ParameterList) );              // ±1 issue
%ignore Tpetra::Import::Import(
        const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& source,
        const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& target,
        const Teuchos::ArrayView<int> & remotePIDs,
        const Teuchos::ArrayView<const LocalOrdinal> & userExportLIDs,
        const Teuchos::ArrayView<const int> & userExportPIDs,
        const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null,
        const Teuchos::RCP<Teuchos::FancyOStream>& out = Teuchos::null); // ±1 issue, needs Teuchos::FancyOStream
%ignore Import (const Teuchos::RCP<const Map<LocalOrdinal, GlobalOrdinal, Node> >& sourceMap,
            const GlobalOrdinal targetMapRemoteOrPermuteGlobalIndices[],
            const int targetMapRemoteOrPermuteProcessRanks[],
            const LocalOrdinal numTargetMapRemoteOrPermuteGlobalIndices,
            const bool mayReorderTargetMapIndicesLocally,
            const Teuchos::RCP<Teuchos::ParameterList>& plist = Teuchos::null,
            const Teuchos::RCP<Teuchos::FancyOStream>& out = Teuchos::null);
%ignore Tpetra::Import::getPermuteFromLIDs;     // ±1 issue
%ignore Tpetra::Import::getPermuteToLIDs;       // ±1 issue
%ignore Tpetra::Import::getRemoteLIDs;          // ±1 issue
%ignore Tpetra::Import::getExportLIDs;          // ±1 issue
%ignore Tpetra::Import::getExportPIDs;          // ±1 issue
%ignore Tpetra::Import::findUnionTargetGIDs;    // ±1 issue
%ignore Tpetra::Import::getDistributor;         // needs Tpetra::Distributor

// Few functions from the base class of Import and Export (Tpetra::Details::Transfer)
%define %tpetra_extend_with_transfer(CLS...)
%extend CLS {
  size_t getNumSameIDs() const {
    return $self->getNumSameIDs();
  }
  size_t getNumPermuteIDs () const {
    return $self->getNumPermuteIDs();
  }
  size_t getNumRemoteIDs () const {
    return $self->getNumRemoteIDs();
  }
  size_t getNumExportIDs () const {
    return $self->getNumExportIDs();
  }
  Teuchos::RCP<const map_type> getSourceMap () const {
    return $self->getSourceMap();
  }
  Teuchos::RCP<const map_type> getTargetMap () const {
    return $self->getTargetMap();
  }
}
%enddef

%tpetra_extend_with_transfer(Tpetra::Import<LO,GO,NO>)
%tpetra_extend_with_transfer(Tpetra::Export<LO,GO,NO>)

// =======================================================================
// Instantiate
// =======================================================================

// First, declare templates
%include "Tpetra_Import_decl.hpp"
%include "Tpetra_Export_decl.hpp"

// Then, label the classes as RCPs
%teuchos_rcp(Tpetra::Export<LO,GO,NO>)
%teuchos_rcp(Tpetra::Import<LO,GO,NO>)

// Finally, instantiate (since they each have the other type as constructors)
%template(TpetraImport) Tpetra::Import<LO,GO,NO>;
%template(TpetraExport) Tpetra::Export<LO,GO,NO>;

// =======================================================================
// Macros
// =======================================================================

%define %tpetra_extend_with_import_export(CLS...)
%extend CLS {
    void doImport (const CLS &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      $self->doImport(source, importer, CM);
    }
    void doImport (const CLS &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      $self->doImport(source, exporter, CM);
    }
    void doExport (const CLS &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      $self->doExport(source, exporter, CM);
    }
    void doExport (const CLS &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      $self->doExport(source, importer, CM);
    }
}
%enddef

