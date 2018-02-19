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
#include "Tpetra_CrsGraph.hpp"
%}

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
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow, \
        const ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::DualView
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const Kokkos::DualView< const size_t *, execution_space > &numEntPerRow, \
        ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::DualView
%ignore Tpetra::CrsGraph::CrsGraph(const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const typename local_graph_type::row_map_type &rowPointers, \
        const typename local_graph_type::entries_type::non_const_type &columnIndices, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::View; ±1 issue
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const local_graph_type &lclGraph, \
        const Teuchos::RCP< Teuchos::ParameterList > &params);                  // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::CrsGraph (const local_graph_type &lclGraph, \
        const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap=Teuchos::null, \
        const Teuchos::RCP< const map_type > &domainMap=Teuchos::null, \
        const Teuchos::RCP< const map_type > &rangeMap=Teuchos::null, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::insertLocalIndices(const LocalOrdinal localRow, const LocalOrdinal numEnt, const LocalOrdinal inds[]);    // prefer ArrayView variant
%ignore Tpetra::CrsGraph::describe;                           // needs Teuchos::FancyOStream
%ignore Tpetra::CrsGraph::getLocalDiagOffsets (const Kokkos::View< size_t *, device_type, Kokkos::MemoryUnmanaged > &offsets) const;    // needs Kokkos::View
%ignore Tpetra::CrsGraph::getNumEntriesPerLocalRowUpperBound; // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsGraph::getLocalGraph;                      // needs Kokkos::StaticCrsGraph
%ignore Tpetra::CrsGraph::getLocalRowView;                    // ±1 issue

// =======================================================================
// Fix ±1 issues
// =======================================================================
%typemap(in)  int localRow %{$1 = *$input - 1;%}
%typemap(argout) const Teuchos::ArrayView< int > &indices %{
  for (int i = 0; i < $1->size(); i++)
    (*$1)[i]++;
%}
%typemap(in, noblock=1) const Teuchos::ArrayView< int > &indices () {

  // Original typemap: convert void* to thinvec reference
  $1 = %static_cast(%const_cast($input, void*), $1_ltype);

  // Construct temporary array and view
  Teuchos::Array<int> tmp = Teuchos::Array<int>($1->size());
  for (int i = 0; i < $1->size(); i++)
    tmp[i] = (*$1)[i] - 1;
  Teuchos::ArrayView<int> tmpview = tmp();

  // Make the input argument point to our temporary vector
  $1 = &tmpview;
}

// =======================================================================
// Make interface more Fortran friendly
// =======================================================================
%extend Tpetra::CrsGraph<LO,GO,NO> {
    CrsGraph(const Teuchos::RCP< const map_type > &rowMap, \
        std::pair<const size_t*,size_t> numEntPerRow, \
        const ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      Teuchos::ArrayRCP<const size_t> numEntPerRowRCP(numEntPerRow.first, 0, numEntPerRow.second, false/*has_ownership*/);
      return new Tpetra::CrsGraph<LO,GO,NO>(rowMap, numEntPerRowRCP, pftype, params);
    }
    CrsGraph(const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        std::pair<const size_t*,size_t> numEntPerRow, \
        const ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      Teuchos::ArrayRCP<const size_t> numEntPerRowRCP(numEntPerRow.first, 0, numEntPerRow.second, false/*has_ownership*/);
      return new Tpetra::CrsGraph<LO,GO,NO>(rowMap, colMap, numEntPerRowRCP, pftype, params);
    }
    CrsGraph(const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        std::pair<size_t*,size_t> rowPointers, \
        std::pair<LO*,size_t> columnIndices, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null) {
      Teuchos::Array<size_t> rowPointersArray(rowPointers.second);
      for (size_t i = 0; i < rowPointers.second; i++)
        rowPointersArray[i] = rowPointers.first[i]-1;
      Teuchos::Array<LO> columnIndicesArray(columnIndices.second);
      for (size_t i = 0; i < columnIndices.second; i++)
        columnIndicesArray[i] = columnIndices.first[i]-1;
      return new Tpetra::CrsGraph<LO,GO,NO>(rowMap, colMap,
        Teuchos::arcpFromArray(rowPointersArray), Teuchos::arcpFromArray(columnIndicesArray), params);
    }
    void insertGlobalIndices(const GO globalRow, std::pair<const GO*,size_t> indices) {
      Teuchos::ArrayView<const GO> indicesView = Teuchos::arrayView(indices.first, indices.second);
      self->insertGlobalIndices(globalRow, indicesView);
    }
    void insertLocalIndices(const LO localRow, std::pair<const LO*,size_t> indices) {
      Teuchos::Array<LO> indicesArray(indices.second);
      for (size_t i = 0; i < indicesArray.size(); i++)
        indicesArray[i] = indices.first[i]-1;
      self->insertLocalIndices(localRow, indicesArray);
    }
    void getGlobalRowCopy(GO GlobalRow, std::pair<GO*,size_t> Indices, size_t &NumIndices) const {
      Teuchos::ArrayView<GO> IndicesView = Teuchos::arrayView(Indices.first, Indices.second);
      self->getGlobalRowCopy(GlobalRow, IndicesView, NumIndices);
    }
    void getLocalRowCopy(LO localRow, std::pair<LO*,size_t> indices, size_t &NumIndices) const {
      Teuchos::ArrayView<LO> indicesView = Teuchos::arrayView(indices.first, indices.second);
      self->getLocalRowCopy(localRow, indicesView, NumIndices);
      for (int i = 0; i < indicesView.size(); i++)
        indicesView[i]++;
    }
    std::pair<const GO*,size_t> getGlobalRowView(const GO gblRow) const {
      Teuchos::ArrayView<const GO> lclColIndsView;
      self->getGlobalRowView(gblRow, lclColIndsView);
      return std::make_pair<const GO*,size_t>(lclColIndsView.getRawPtr(), lclColIndsView.size());
    }
    void setAllIndices(std::pair<size_t*,size_t> rowPointers, std::pair<LO*,size_t> columnIndices, std::pair<SC*,size_t> val) {
      Teuchos::ArrayRCP<size_t> rowPointersArrayRCP(rowPointers.second);
      for (int i = 0; i < rowPointersArrayRCP.size(); i++)
        rowPointersArrayRCP[i] = rowPointers.first[i]-1;
      Teuchos::ArrayRCP<LO> columnIndicesArrayRCP(columnIndices.second);
      for (int i = 0; i < columnIndicesArrayRCP.size(); i++)
        columnIndicesArrayRCP[i] = columnIndices.first[i]-1;
      self->setAllIndices(rowPointersArrayRCP, columnIndicesArrayRCP);
    }
    // NOTE: This is semantically different function from Tpetra. Here, we *require* that user already allocated the arrays to store the data
    void getNodeRowPtrs(std::pair<size_t*,size_t> rowPointers) const {
      auto rowPointersArrayRCP = self->getNodeRowPtrs();
      TEUCHOS_TEST_FOR_EXCEPTION(rowPointersArrayRCP.size() != rowPointers.second, std::runtime_error, "Wrong rowPointers size");
      auto n = rowPointers.second;
      for (int i = 0; i < n; i++)
        rowPointers.first[i] = rowPointersArrayRCP[i]+1;
    }
    void getNodePackedIndices(std::pair<size_t*,size_t> columnIndices) const {
      auto columnIndicesArrayRCP = self->getNodeRowPtrs();
      TEUCHOS_TEST_FOR_EXCEPTION(columnIndicesArrayRCP.size() != columnIndices.second, std::runtime_error, "Wrong columnIndices size");
      auto nnz = columnIndices.second;
      for (int i = 0; i < nnz; i++)
        columnIndices.first[i] = columnIndicesArrayRCP[i]+1;
    }
    void getLocalDiagOffsets (std::pair<size_t*,size_t> offsets) const {
      TEUCHOS_TEST_FOR_EXCEPTION(self->getNodeNumRows() != offsets.second, std::runtime_error, "Wrong offsets size");
      Teuchos::ArrayRCP<size_t> offsetsArrayRCP(offsets.first, 0, offsets.second, false/*has_ownership*/);
      self->getLocalDiagOffsets(offsetsArrayRCP);
    }
    void doImport (const Tpetra::CrsGraph<LO,GO,NO> &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      self->doImport(source, importer, CM);
    }
    void doImport (const Tpetra::CrsGraph<LO,GO,NO> &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      self->doImport(source, exporter, CM);
    }
    void doExport (const Tpetra::CrsGraph<LO,GO,NO> &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      self->doExport(source, exporter, CM);
    }
    void doExport (const Tpetra::CrsGraph<LO,GO,NO> &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      self->doExport(source, importer, CM);
    }
}
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::ArrayRCP< const size_t > &numEntPerRow, \
        const ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const Teuchos::ArrayRCP< const size_t > &numEntPerRow, \
        ProfileType pftype=DynamicProfile, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);
%ignore Tpetra::CrsGraph::CrsGraph (const Teuchos::RCP< const map_type > &rowMap, \
        const Teuchos::RCP< const map_type > &colMap, \
        const Teuchos::ArrayRCP< size_t > &rowPointers, \
        const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices, \
        const Teuchos::RCP< Teuchos::ParameterList > &params=Teuchos::null);    // needs Teuchos::ArrayRCP; ±1 issue
%ignore Tpetra::CrsGraph::insertGlobalIndices (const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &indices);
%ignore Tpetra::CrsGraph::insertLocalIndices(const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &indices);
%ignore Tpetra::CrsGraph::getGlobalRowCopy(GlobalOrdinal GlobalRow, const Teuchos::ArrayView< GlobalOrdinal > &Indices, size_t &NumIndices) const;
%ignore Tpetra::CrsGraph::getLocalRowCopy(LocalOrdinal LocalRow, const Teuchos::ArrayView< LocalOrdinal > &indices, size_t &NumIndices) const;
%ignore Tpetra::CrsGraph::getGlobalRowView(const GlobalOrdinal gblRow, Teuchos::ArrayView< const GlobalOrdinal > &gblColInds) const;
%ignore Tpetra::CrsGraph::getLocalRowView(const LocalOrdinal lclRow, Teuchos::ArrayView< const LocalOrdinal > &lclColInds) const;
%ignore Tpetra::CrsGraph::getNodeRowPtrs() const;
%ignore Tpetra::CrsGraph::getNodePackedIndices() const;
%ignore Tpetra::CrsGraph::getLocalDiagOffsets (Teuchos::ArrayRCP< size_t > &offsets) const;
%ignore Tpetra::CrsGraph::setAllIndices (const Teuchos::ArrayRCP< size_t > &rowPointers, const Teuchos::ArrayRCP< LocalOrdinal > &columnIndices);
%ignore Tpetra::CrsGraph::setAllIndices (const typename local_graph_type::row_map_type &rowPointers, const typename local_graph_type::entries_type::non_const_type &columnIndices);


%teuchos_rcp(Tpetra::CrsGraph<LO,GO,NO,NO::classic>)

%include "Tpetra_CrsGraph_decl.hpp"

%template(TpetraCrsGraph) Tpetra::CrsGraph<LO,GO,NO>;
