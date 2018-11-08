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
#include "Tpetra_CrsMatrix.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================
// Ignore Details namespace
%ignore Details;
%ignore Tpetra::CrsMatrix::pack_functor;
%ignore Tpetra::CrsMatrix::getNode;
// Ignore implementation details
%ignore Tpetra::CrsMatrix::checkSizes;
%ignore Tpetra::CrsMatrix::copyAndPermute;
%ignore Tpetra::CrsMatrix::copyAndPermuteNew;
%ignore Tpetra::CrsMatrix::pack;
%ignore Tpetra::CrsMatrix::packAndPrepare;
%ignore Tpetra::CrsMatrix::packAndPrepareNew;
%ignore Tpetra::CrsMatrix::packNew;
%ignore Tpetra::CrsMatrix::packNonStatic;
%ignore Tpetra::CrsMatrix::unpackAndCombine;
%ignore Tpetra::CrsMatrix::unpackAndCombineNew;
%ignore Tpetra::CrsMatrix::useNewInterface;

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const typename local_matrix_type::row_map_type& rowPointers,
               const typename local_graph_type::entries_type::non_const_type& columnIndices,
               const typename local_matrix_type::values_type& values,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs Kokkos containers
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const local_matrix_type& lclMatrix,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs Kokkos containers
%ignore Tpetra::CrsMatrix::CrsMatrix (const local_matrix_type& lclMatrix,
               const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap = Teuchos::null,
               const Teuchos::RCP<const map_type>& domainMap = Teuchos::null,
               const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs Kokkos containers
%ignore Tpetra::CrsMatrix::setAllValues(const typename local_matrix_type::row_map_type& ptr,
               const typename local_graph_type::entries_type::non_const_type& ind,
               const typename local_matrix_type::values_type& val);                         // needs Kokkos::View
%ignore Tpetra::CrsMatrix::insertGlobalValues (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const Scalar vals[], const GlobalOrdinal inds[]); // prefer ArrayView version
%ignore Tpetra::CrsMatrix::insertLocalValues (const LocalOrdinal localRow, const LocalOrdinal numEnt, const Scalar vals[], const LocalOrdinal cols[]); // prefer ArrayView version
%ignore Tpetra::CrsMatrix::replaceGlobalValues (const GlobalOrdinal globalRow, const typename UnmanagedView< GlobalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals) const;  // needs Kokkos::UnmanagedView
%ignore Tpetra::CrsMatrix::replaceGlobalValues (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const Scalar vals[], const GlobalOrdinal cols[]) const;                                                         // prefer ArrayView variant
%ignore Tpetra::CrsMatrix::sumIntoGlobalValues (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const Scalar vals[], const GlobalOrdinal cols[], const bool atomic=useAtomicUpdatesByDefault); // prefer ArrayView variant
%ignore Tpetra::CrsMatrix::sumIntoLocalValues (const LocalOrdinal localRow, const typename UnmanagedView< LocalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals, const bool atomic=useAtomicUpdatesByDefault) const; // needs Kokkos::UnmanagedView
%ignore Tpetra::CrsMatrix::sumIntoLocalValues (const LocalOrdinal localRow, const LocalOrdinal numEnt, const Scalar vals[], const LocalOrdinal cols[], const bool atomic=useAtomicUpdatesByDefault) const; // prefer ArrayView variant
%ignore Tpetra::CrsMatrix::operator();                      // needs operator() treatment
%ignore Tpetra::CrsMatrix::describe;                        // needs Teuchos::FancyOStream
%ignore Tpetra::CrsMatrix::add;                             // needs Tpetra::RowMatrix
%ignore Tpetra::CrsMatrix::getGraph;                        // needs Tpetra::RowGraph
%ignore Tpetra::CrsMatrix::getLocalDiagCopy;                // needs Tpetra::Vector
%ignore Tpetra::CrsMatrix::getLocalMatrix;                  // needs KokkosSparse::CrsMatrix
%ignore Tpetra::CrsMatrix::getLocalValuesView;              // needs Kokkos::View
%ignore Tpetra::CrsMatrix::leftScale;                       // needs Tpetra::Vector
%ignore Tpetra::CrsMatrix::rightScale;                      // needs Tpetra::Vector
%ignore Tpetra::CrsMatrix::reorderedGaussSeidel;
%ignore Tpetra::CrsMatrix::reorderedGaussSeidelCopy;

// =======================================================================
// Fix ±1 issues
// =======================================================================
%typemap(in)  int localRow %{$1 = *$input - 1;%}            // int = LO
%ignore Tpetra::CrsMatrix::getLocalRowView;                 // ±1 issue
%ignore Tpetra::CrsMatrix::reorderedLocalGaussSeidel;       // ±1 issue
%ignore Tpetra::CrsMatrix::replaceLocalValues;              // ±1 issue
%ignore Tpetra::CrsMatrix::sumIntoLocalValues;              // ±1 issue
%ignore Tpetra::CrsMatrix::transformLocalValues;            // ±1 issue

// =======================================================================
// Make interface more Fortran friendly
// =======================================================================
%extend Tpetra::CrsMatrix<SC,LO,GO,NO> {
    CrsMatrix(const Teuchos::RCP<const map_type>& rowMap, std::pair<const size_t*,size_t> NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile, const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      Teuchos::ArrayRCP<const size_t> NumEntriesPerRowToAllocArrayRCP(NumEntriesPerRowToAlloc.first, 0, NumEntriesPerRowToAlloc.second, false/*has_ownership*/);
      return new Tpetra::CrsMatrix<SC,LO,GO,NO>(rowMap, NumEntriesPerRowToAllocArrayRCP, pftype, params);
    }
    CrsMatrix(const Teuchos::RCP<const map_type>& rowMap, const Teuchos::RCP<const map_type>& colMap, std::pair<const size_t*,size_t> NumEntriesPerRowToAlloc, ProfileType pftype = DynamicProfile, const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      Teuchos::ArrayRCP<const size_t> NumEntriesPerRowToAllocArrayRCP(NumEntriesPerRowToAlloc.first, 0, NumEntriesPerRowToAlloc.second, false/*has_ownership*/);
      return new Tpetra::CrsMatrix<SC,LO,GO,NO>(rowMap, colMap, NumEntriesPerRowToAllocArrayRCP, pftype, params);
    }
    CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, const Teuchos::RCP<const map_type>& colMap, std::pair<size_t*,size_t> rowPointers, std::pair<LO*,size_t> columnIndices, std::pair<SC*,size_t> values, const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null) {
      Teuchos::ArrayRCP<size_t> rowPointersArrayRCP(rowPointers.second);
      for (int i = 0; i < rowPointersArrayRCP.size(); i++)
        rowPointersArrayRCP[i] = rowPointers.first[i]-1;
      Teuchos::ArrayRCP<LO> columnIndicesArrayRCP(columnIndices.second);
      for (int i = 0; i < columnIndicesArrayRCP.size(); i++)
        columnIndicesArrayRCP[i] = columnIndices.first[i]-1;
      Teuchos::ArrayRCP<SC> valuesArrayRCP(values.first, 0, values.second, false/*has_ownership*/);
      return new Tpetra::CrsMatrix<SC,LO,GO,NO>(rowMap, colMap, rowPointersArrayRCP, columnIndicesArrayRCP, valuesArrayRCP, params);
    }
    void insertGlobalValues(const GO globalRow, std::pair<const GO*,size_t> cols, std::pair<const SC*,size_t> vals) {
      Teuchos::ArrayView<const GO> colsView = Teuchos::arrayView(cols.first, cols.second);
      Teuchos::ArrayView<const SC> valsView = Teuchos::arrayView(vals.first, vals.second);
      self->insertGlobalValues(globalRow, colsView, valsView);
    }
    void insertLocalValues(const LO localRow, std::pair<const LO*,size_t> cols, std::pair<const SC*,size_t> vals) {
      Teuchos::Array<LO> colsArray(cols.second);
      for (int i = 0; i < colsArray.size(); i++)
        colsArray[i] = cols.first[i] - 1;
      Teuchos::ArrayView<const SC> valsView = Teuchos::arrayView(vals.first, vals.second);
      self->insertLocalValues(localRow, colsArray, valsView);
    }
    LO replaceGlobalValues(const GO globalRow, std::pair<const GO*,size_t> cols, std::pair<const SC*,size_t> vals) const {
      Teuchos::ArrayView<const GO> colsView = Teuchos::arrayView(cols.first, cols.second);
      Teuchos::ArrayView<const SC> valsView = Teuchos::arrayView(vals.first, vals.second);
      return self->replaceGlobalValues(globalRow, colsView, valsView);
    }
    LO replaceLocalValues(const LO localRow, std::pair<const LO*,size_t> cols, std::pair<const SC*,size_t> vals) const {
      Teuchos::Array<LO> colsArray(cols.second);
      for (int i = 0; i < colsArray.size(); i++)
        colsArray[i] = cols.first[i] - 1;
      Teuchos::ArrayView<const SC> valsView = Teuchos::arrayView(vals.first, vals.second);
      return self->replaceLocalValues(localRow, colsArray, valsView);
    }
    LO sumIntoGlobalValues(const GO globalRow, std::pair<const GO*,size_t> cols, std::pair<const SC*,size_t> vals) {
      Teuchos::ArrayView<const GO> colsView = Teuchos::arrayView(cols.first, cols.second);
      Teuchos::ArrayView<const SC> valsView = Teuchos::arrayView(vals.first, vals.second);
      return self->sumIntoGlobalValues(globalRow, colsView, valsView, false); // TODO: for now, we only run in serial, no atomics necessary
    }
    LO sumIntoLocalValues(const LO localRow, std::pair<const LO*,size_t> cols, std::pair<const SC*,size_t> vals) const {
      Teuchos::Array<LO> colsArray(cols.second);
      for (int i = 0; i < colsArray.size(); i++)
        colsArray[i] = cols.first[i] - 1;
      Teuchos::ArrayView<const SC> valsView = Teuchos::arrayView(vals.first, vals.second);
      return self->sumIntoLocalValues(localRow, colsArray, valsView, false/*atomic*/); // TODO: for now, we only run in serial, no atomics necessary
    }
    void setAllValues(std::pair<size_t*,size_t> ptr, std::pair<LO*,size_t> ind, std::pair<SC*,size_t> val) {
      Teuchos::ArrayRCP<size_t> ptrArrayRCP(ptr.second);
      for (int i = 0; i < ptrArrayRCP.size(); i++)
        ptrArrayRCP[i] = ptr.first[i]-1;
      Teuchos::ArrayRCP<LO> indArrayRCP(ind.second);
      for (int i = 0; i < indArrayRCP.size(); i++)
        indArrayRCP[i] = ind.first[i]-1;
      Teuchos::ArrayRCP<SC> valArrayRCP(val.first, 0, val.second, false/*has_ownership*/);
      self->setAllValues(ptrArrayRCP, indArrayRCP, valArrayRCP);
    }
    // NOTE: This is semantically different function from Tpetra. Here, we *require* that user already allocated the arrays to store the data
    void getAllValues(std::pair<size_t*,size_t> rowPointers, std::pair<LO*,size_t> columnIndices, std::pair<SC*,size_t> values) const {
      Teuchos::ArrayRCP<const size_t> rowPointersArrayRCP;
      Teuchos::ArrayRCP<const LO>     columnIndicesArrayRCP;
      Teuchos::ArrayRCP<const SC>     valuesArrayRCP;
      self->getAllValues(rowPointersArrayRCP, columnIndicesArrayRCP, valuesArrayRCP);
      TEUCHOS_TEST_FOR_EXCEPTION(rowPointersArrayRCP.size()   != rowPointers.second,    std::runtime_error, "Wrong rowPointers size");
      TEUCHOS_TEST_FOR_EXCEPTION(columnIndicesArrayRCP.size() != columnIndices.second,  std::runtime_error, "Wrong columnIndices size");
      TEUCHOS_TEST_FOR_EXCEPTION(valuesArrayRCP.size()        != values.second,         std::runtime_error, "Wrong values size");
      auto n = rowPointers.second;
      for (int i = 0; i < n; i++)
        rowPointers.first[i] = rowPointersArrayRCP[i]+1;
      auto nnz = columnIndices.second;
      for (int i = 0; i < nnz; i++) {
        columnIndices.first[i] = columnIndicesArrayRCP[i]+1;
        values       .first[i] = valuesArrayRCP[i];
      }
    }
    void getGlobalRowCopy(GO GlobalRow, std::pair<GO*,size_t> Indices, std::pair<SC*, size_t> Values, size_t &NumIndices) const {
      Teuchos::ArrayView<GO> IndicesView = Teuchos::arrayView(Indices.first, Indices.second);
      Teuchos::ArrayView<SC> ValuesView  = Teuchos::arrayView(Values.first, Values.second);

      self->getGlobalRowCopy(GlobalRow, IndicesView, ValuesView, NumIndices);
    }
    void getLocalRowCopy(LO localRow, std::pair<LO*,size_t> colInds, std::pair<SC*, size_t> vals, size_t &NumIndices) const {
      Teuchos::ArrayView<LO> colIndsView = Teuchos::arrayView(colInds.first, colInds.second);
      Teuchos::ArrayView<SC> valsView  = Teuchos::arrayView(vals.first, vals.second);

      self->getLocalRowCopy(localRow, colIndsView, valsView, NumIndices);

      for (int i = 0; i < colIndsView.size(); i++)
        colIndsView[i]++;
    }
    void getGlobalRowView(GO GlobalRow, std::pair<const GO*,size_t> Indices, std::pair<const SC*, size_t> values) const {
      Teuchos::ArrayView<const GO> IndicesView = Teuchos::arrayView(Indices.first, Indices.second);
      Teuchos::ArrayView<const SC> valuesView  = Teuchos::arrayView(values.first, values.second);

      self->getGlobalRowView(GlobalRow, IndicesView, valuesView);
    }
    void getLocalDiagOffsets (std::pair<size_t*,size_t> offsets) const {
      TEUCHOS_TEST_FOR_EXCEPTION(self->getNodeNumRows() != offsets.second, std::runtime_error, "Wrong offsets size");
      Teuchos::ArrayRCP<size_t> offsetsArrayRCP(offsets.first, 0, offsets.second, false/*has_ownership*/);
      self->getLocalDiagOffsets(offsetsArrayRCP);
    }
    void doImport (const Tpetra::CrsMatrix<SC,LO,GO,NO> &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      self->doImport(source, importer, CM);
    }
    void doImport (const Tpetra::CrsMatrix<SC,LO,GO,NO> &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      self->doImport(source, exporter, CM);
    }
    void doExport (const Tpetra::CrsMatrix<SC,LO,GO,NO> &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      self->doExport(source, exporter, CM);
    }
    void doExport (const Tpetra::CrsMatrix<SC,LO,GO,NO> &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      self->doExport(source, importer, CM);
    }
}
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs ArrayRCP
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc,
               ProfileType pftype = DynamicProfile,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs ArrayRCP
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap,
               const Teuchos::RCP<const map_type>& colMap,
               const Teuchos::ArrayRCP<size_t>& rowPointers,
               const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices,
               const Teuchos::ArrayRCP<Scalar>& values,
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsMatrix::insertGlobalValues(const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals);
%ignore Tpetra::CrsMatrix::insertLocalValues (const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals);
%ignore Tpetra::CrsMatrix::getGlobalRowCopy(GlobalOrdinal GlobalRow, const Teuchos::ArrayView< GlobalOrdinal > &Indices, const Teuchos::ArrayView< Scalar > &Values, size_t &NumEntries) const;
%ignore Tpetra::CrsMatrix::replaceGlobalValues(const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals) const;
%ignore Tpetra::CrsMatrix::replaceLocalValues(const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals) const;
%ignore Tpetra::CrsMatrix::sumIntoGlobalValues(const GlobalOrdinal globalRow, const Teuchos::ArrayView< const GlobalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals, const bool atomic=useAtomicUpdatesByDefault);
%ignore Tpetra::CrsMatrix::setAllValues(const Teuchos::ArrayRCP< size_t > &ptr, const Teuchos::ArrayRCP< LocalOrdinal > &ind, const Teuchos::ArrayRCP< Scalar > &val);
%ignore Tpetra::CrsMatrix::getAllValues(Teuchos::ArrayRCP< const size_t > &rowPointers, Teuchos::ArrayRCP< const LocalOrdinal > &columnIndices, Teuchos::ArrayRCP< const Scalar > &values) const;
%ignore Tpetra::CrsMatrix::getGlobalRowCopy(GlobalOrdinal GlobalRow, const Teuchos::ArrayView< GlobalOrdinal > &Indices, const Teuchos::ArrayView< Scalar > &Values, size_t &NumEntries) const;
%ignore Tpetra::CrsMatrix::getLocalRowCopy(LocalOrdinal localRow, const Teuchos::ArrayView< LocalOrdinal > &colInds, const Teuchos::ArrayView< Scalar > &vals, size_t &numEntries) const;
%ignore Tpetra::CrsMatrix::getGlobalRowView(GlobalOrdinal GlobalRow, Teuchos::ArrayView< const GlobalOrdinal > &indices, Teuchos::ArrayView< const Scalar > &values) const;
%ignore Tpetra::CrsMatrix::getLocalRowView(LocalOrdinal LocalRow, Teuchos::ArrayView< const LocalOrdinal > &indices, Teuchos::ArrayView< const Scalar > &values) const;
%ignore Tpetra::CrsMatrix::getLocalRowViewRaw(const LocalOrdinal lclRow, LocalOrdinal &numEnt, const LocalOrdinal *&lclColInds, const Scalar *&vals) const;
%ignore Tpetra::CrsMatrix::getAllValues (Teuchos::ArrayRCP<const size_t>& rowPointers, Teuchos::ArrayRCP<const LocalOrdinal>& columnIndices, Teuchos::ArrayRCP<const Scalar>& values) const;
%ignore Tpetra::CrsMatrix::getLocalDiagOffsets;             // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsMatrix::sumIntoLocalValues (const LocalOrdinal localRow, const typename UnmanagedView< LocalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals, const bool atomic=useAtomicUpdatesByDefault) const;
%ignore Tpetra::CrsMatrix::sumIntoLocalValues (const LocalOrdinal localRow, const Teuchos::ArrayView< const LocalOrdinal > &cols, const Teuchos::ArrayView< const Scalar > &vals, const bool atomic=useAtomicUpdatesByDefault) const;
%ignore Tpetra::CrsMatrix::sumIntoLocalValues (const LocalOrdinal localRow, const LocalOrdinal numEnt, const Scalar vals[], const LocalOrdinal cols[], const bool atomic=useAtomicUpdatesByDefault);


%include "Tpetra_CrsMatrix_decl.hpp"

%teuchos_rcp(Tpetra::CrsMatrix<SC,LO,GO,NO>)
%template(TpetraCrsMatrix) Tpetra::CrsMatrix<SC,LO,GO,NO>;

// Operator to Matrix conversion
%{
#include "Tpetra_CrsMatrix.hpp"
%}
%inline %{
namespace ForTrilinos {
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  operator_to_matrix(Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > op) {
    auto A = Teuchos::rcp_dynamic_cast<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(op);
    TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), std::runtime_error, "operator_to_matrix: the provided operator is not a Tpetra CrsMatrix");
    return A;
  }
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  matrix_to_operator(Teuchos::RCP<Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> > A) {
    auto op = Teuchos::rcp_dynamic_cast<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> >(A);
    return op;
  }
}
%}
%template(operator_to_matrix) ForTrilinos::operator_to_matrix<SC,LO,GO,NO>;
%template(matrix_to_operator) ForTrilinos::matrix_to_operator<SC,LO,GO,NO>;
