/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"
%import <Teuchos_Comm.i>

%{
#include "Tpetra_CrsGraph.hpp"
%}

%ignore Tpetra::Details::HasDeprecatedMethods2630_WarningThisClassIsNotForUsers;

// =======================================================================
// Ignore permanently
// =======================================================================
%ignore Tpetra::CrsGraph::getNode;
%ignore Tpetra::CrsGraph::checkSizes;
%ignore Tpetra::CrsGraph::copyAndPermute;
%ignore Tpetra::CrsGraph::packAndPrepare;
%ignore Tpetra::CrsGraph::pack;
%ignore Tpetra::CrsGraph::pack_functor;
%ignore Tpetra::CrsGraph::unpackAndCombine;
%ignore Tpetra::CrsGraph::SLocalGlobalViews;
%ignore Tpetra::CrsGraph::SLocalGlobalNCViews;

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap,
        const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow,
        const ProfileType pftype=DynamicProfile,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::DualView
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP<const map_type>& rowMap,
        const Teuchos::RCP<const map_type>& colMap,
        const Kokkos::DualView<const size_t*, execution_space>& numEntPerRow,
        const ProfileType pftype = DynamicProfile,
        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);    // needs Kokkos::DualView
%ignore Tpetra::CrsGraph::CrsGraph(const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::RCP< const map_type > &colMap,
        const typename local_graph_type::row_map_type &rowPointers,
        const typename local_graph_type::entries_type::non_const_type &columnIndices,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::View; ±1 issue
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::RCP< const map_type > &colMap,
        const local_graph_type &lclGraph,
        const Teuchos::RCP< Teuchos::ParameterList > &params);                  // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::CrsGraph (const local_graph_type &lclGraph,
        const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::RCP< const map_type > &colMap=Teuchos::null,
        const Teuchos::RCP< const map_type > &domainMap=Teuchos::null,
        const Teuchos::RCP< const map_type > &rangeMap=Teuchos::null,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::insertLocalIndices(const LocalOrdinal localRow, const LocalOrdinal numEnt, const LocalOrdinal inds[]);    // prefer ArrayView variant
%ignore Tpetra::CrsGraph::describe;                           // needs Teuchos::FancyOStream
%ignore Tpetra::CrsGraph::getLocalDiagOffsets (const Kokkos::View< size_t *, device_type, Kokkos::MemoryUnmanaged > &offsets) const;    // needs Kokkos::View
%ignore Tpetra::CrsGraph::getNumEntriesPerLocalRowUpperBound; // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsGraph::getLocalGraph;                      // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::getLocalRowView;                    // ±1 issue

// Returns ArrayView by reference
%ignore Tpetra::CrsGraph::getGlobalRowView;

// =======================================================================
// Fix ±1 issues
// =======================================================================
%apply int INDEX { int localRow };

%apply int { size_t getNumEntriesInLocalRow,
             size_t getNumAllocatedEntriesInLocalRow}

// =======================================================================
// Make interface more Fortran friendly
// =======================================================================
%extend Tpetra::CrsGraph<LO,GO,NO> {
    CrsGraph(const Teuchos::RCP< const map_type > &rowMap,
        Teuchos::ArrayView<const size_t> numEntPerRow,
        const ProfileType pftype=DynamicProfile,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      return new Tpetra::CrsGraph<LO,GO,NO>(rowMap, arcpFromArrayView(numEntPerRow), pftype, params);
    }
    CrsGraph(const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::RCP< const map_type > &colMap,
        Teuchos::ArrayView<const size_t> numEntPerRow,
        const ProfileType pftype=DynamicProfile,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      return new Tpetra::CrsGraph<LO,GO,NO>(rowMap, colMap, arcpFromArrayView(numEntPerRow), pftype, params);
    }
    CrsGraph(const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::RCP< const map_type > &colMap,
        Teuchos::ArrayView<size_t> rowPointers,
        Teuchos::ArrayView<LO> columnIndices,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      Teuchos::Array<size_t> rowPointersArray(rowPointers.size());
      for (size_t i = 0; i < rowPointers.size(); i++)
        rowPointersArray[i] = rowPointers[i]-1;
      Teuchos::Array<LO> columnIndicesArray(columnIndices.size());
      for (size_t i = 0; i < columnIndices.size(); i++)
        columnIndicesArray[i] = columnIndices[i]-1;
      return new Tpetra::CrsGraph<LO,GO,NO>(rowMap, colMap,
        Teuchos::arcpFromArray(rowPointersArray), Teuchos::arcpFromArray(columnIndicesArray), params);
    }
    void insertGlobalIndices(const GO globalRow, Teuchos::ArrayView<const GO> indicesView) {
      $self->insertGlobalIndices(globalRow, indicesView);
    }
    void insertLocalIndices(const LO localRow, Teuchos::ArrayView<const LO> indices) {
      Teuchos::Array<LO> indicesArray(indices.size());
      for (size_t i = 0; i < indicesArray.size(); i++)
        indicesArray[i] = indices[i]-1;
      $self->insertLocalIndices(localRow, indicesArray);
    }
    void getLocalRowCopy(LO localRow, Teuchos::ArrayView<LO> indicesView, size_t &NumIndices) const {
      $self->getLocalRowCopy(localRow, indicesView, NumIndices);
      for (int i = 0; i < indicesView.size(); i++)
        indicesView[i]++;
    }
    void setAllIndices(Teuchos::ArrayView<size_t> rowPointers, Teuchos::ArrayView<LO> columnIndices, Teuchos::ArrayView<SC> val) {
      Teuchos::ArrayRCP<size_t> rowPointersArrayRCP(rowPointers.size());
      for (int i = 0; i < rowPointersArrayRCP.size(); i++)
        rowPointersArrayRCP[i] = rowPointers[i]-1;
      Teuchos::ArrayRCP<LO> columnIndicesArrayRCP(columnIndices.size());
      for (int i = 0; i < columnIndicesArrayRCP.size(); i++)
        columnIndicesArrayRCP[i] = columnIndices[i]-1;
      $self->setAllIndices(rowPointersArrayRCP, columnIndicesArrayRCP);
    }
    // NOTE: This is semantically different function from Tpetra. Here, we *require* that user already allocated the arrays to store the data
    void getNodeRowPtrs(Teuchos::ArrayView<size_t> rowPointers) const {
      auto rowPointersArrayRCP = $self->getNodeRowPtrs();
      TEUCHOS_TEST_FOR_EXCEPTION(rowPointersArrayRCP.size() != rowPointers.size(), std::runtime_error, "Wrong rowPointers size");
      auto n = rowPointers.size();
      for (int i = 0; i < n; i++)
        rowPointers[i] = rowPointersArrayRCP[i]+1;
    }
    void getNodePackedIndices(Teuchos::ArrayView<size_t> columnIndices) const {
      auto columnIndicesArrayRCP = $self->getNodeRowPtrs();
      TEUCHOS_TEST_FOR_EXCEPTION(columnIndicesArrayRCP.size() != columnIndices.size(), std::runtime_error, "Wrong columnIndices size");
      auto nnz = columnIndices.size();
      for (int i = 0; i < nnz; i++)
        columnIndices[i] = columnIndicesArrayRCP[i]+1;
    }
    void doImport (const Tpetra::CrsGraph<LO,GO,NO> &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      $self->doImport(source, importer, CM);
    }
    void doImport (const Tpetra::CrsGraph<LO,GO,NO> &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      $self->doImport(source, exporter, CM);
    }
    void doExport (const Tpetra::CrsGraph<LO,GO,NO> &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      $self->doExport(source, exporter, CM);
    }
    void doExport (const Tpetra::CrsGraph<LO,GO,NO> &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      $self->doExport(source, importer, CM);
    }
}
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::ArrayRCP< const size_t > &numEntPerRow,
        const ProfileType pftype=DynamicProfile,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP<const map_type>& rowMap,                                                                                    const Teuchos::RCP<const map_type>& colMap,
        const Teuchos::ArrayRCP<const size_t>& numEntPerRow,
        const ProfileType pftype = DynamicProfile,
        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap,
        const Teuchos::RCP< const map_type > &colMap,
        const Teuchos::ArrayRCP< size_t > &rowPointers,
        const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices,
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Teuchos::ArrayRCP; ±1 issue
%ignore Tpetra::CrsGraph::insertGlobalIndices (const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &indices);
%ignore Tpetra::CrsGraph::insertLocalIndices(const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &indices);
%ignore Tpetra::CrsGraph::getLocalRowCopy(LocalOrdinal LocalRow, const Teuchos::ArrayView< LocalOrdinal > &indices, size_t &NumIndices) const;
%ignore Tpetra::CrsGraph::getNodeRowPtrs() const;
%ignore Tpetra::CrsGraph::getNodePackedIndices() const;
%ignore Tpetra::CrsGraph::getLocalDiagOffsets;
%ignore Tpetra::CrsGraph::setAllIndices (const Teuchos::ArrayRCP< size_t > &rowPointers, const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices);
%ignore Tpetra::CrsGraph::setAllIndices (const typename local_graph_type::row_map_type &rowPointers, const typename local_graph_type::entries_type::non_const_type &columnIndices);


%include "Tpetra_CrsGraph_decl.hpp"

%teuchos_rcp(Tpetra::CrsGraph<LO,GO,NO>)
%template(TpetraCrsGraph) Tpetra::CrsGraph<LO,GO,NO>;
