/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"
%import <Teuchos_ArrayView.i>
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
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap, \
               const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc, \
               ProfileType pftype = DynamicProfile, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs ArrayRCP
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc, \
               ProfileType pftype = DynamicProfile, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs ArrayRCP
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap, \
               const typename local_matrix_type::row_map_type& rowPointers, \
               const typename local_graph_type::entries_type::non_const_type& columnIndices, \
               const typename local_matrix_type::values_type& values, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs Kokkos containers
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap, \
               const Teuchos::ArrayRCP<size_t>& rowPointers, \
               const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices, \
               const Teuchos::ArrayRCP<Scalar>& values, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap, \
               const local_matrix_type& lclMatrix, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs Kokkos containers
%ignore Tpetra::CrsMatrix::CrsMatrix (const local_matrix_type& lclMatrix, \
               const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap = Teuchos::null, \
               const Teuchos::RCP<const map_type>& domainMap = Teuchos::null, \
               const Teuchos::RCP<const map_type>& rangeMap = Teuchos::null, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);         // needs Kokkos containers
%ignore Tpetra::CrsMatrix::setAllValues (const Teuchos::ArrayRCP<size_t>& ptr, \
               const Teuchos::ArrayRCP<LocalOrdinal>& ind, \
               const Teuchos::ArrayRCP<Scalar>& val);                                       // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsMatrix::setAllValues(const typename local_matrix_type::row_map_type& ptr, \
               const typename local_graph_type::entries_type::non_const_type& ind, \
               const typename local_matrix_type::values_type& val);                         // needs Kokkos::View
%ignore Tpetra::CrsMatrix::getAllValues (Teuchos::ArrayRCP<const size_t>& rowPointers, \
               Teuchos::ArrayRCP<const LocalOrdinal>& columnIndices, \
               Teuchos::ArrayRCP<const Scalar>& values) const;                              // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsMatrix::getAllValues (Teuchos::ArrayRCP<const size_t>& rowPointers, \
               Teuchos::ArrayRCP<const LocalOrdinal>& columnIndices, \
               Teuchos::ArrayRCP<const Scalar>& values) const;                              // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsMatrix::operator();                      // needs operator() treatment
%ignore Tpetra::CrsMatrix::describe;                        // needs Teuchos::FancyOStream
%ignore Tpetra::CrsMatrix::add;                             // needs Tpetra::RowMatrix
%ignore Tpetra::CrsMatrix::getGraph;                        // needs Tpetra::RowGraph
%ignore Tpetra::CrsMatrix::getLocalDiagCopy;                // needs Tpetra::Vector
%ignore Tpetra::CrsMatrix::getLocalDiagOffsets;             // needs Teuchos::ArrayRCP
%ignore Tpetra::CrsMatrix::getLocalMatrix;                  // needs KokkosSparse::CrsMatrix
%ignore Tpetra::CrsMatrix::getLocalValuesView;              // needs Kokkos::View
%ignore Tpetra::CrsMatrix::leftScale;                       // needs Tpetra::Vector
%ignore Tpetra::CrsMatrix::rightScale;                      // needs Tpetra::Vector

// =======================================================================
// Fix ±1 issues
// =======================================================================
%typemap(in)  int localRow %{$1 = *$input - 1;%}
%ignore Tpetra::CrsMatrix::getLocalRowCopy;                 // ±1 issue
%ignore Tpetra::CrsMatrix::getLocalRowView;                 // ±1 issue
%ignore Tpetra::CrsMatrix::getNumEntriesInLocalRow;         // ±1 issue
%ignore Tpetra::CrsMatrix::insertLocalValues;               // ±1 issue
%ignore Tpetra::CrsMatrix::reorderedLocalGaussSeidel;       // ±1 issue
%ignore Tpetra::CrsMatrix::replaceLocalValues;              // ±1 issue
%ignore Tpetra::CrsMatrix::sumIntoLocalValues;              // ±1 issue
%ignore Tpetra::CrsMatrix::transformLocalValues;            // ±1 issue

// =======================================================================
// Make interface more Fortran friendly
// =======================================================================
%ignore Tpetra::CrsMatrix::getGlobalRowCopy(GlobalOrdinal GlobalRow, const Teuchos::ArrayView< GlobalOrdinal > &Indices, const Teuchos::ArrayView< Scalar > &Values, size_t &NumEntries) const;
%extend Tpetra::CrsMatrix<SC,LO,GO,NO,false> {
    void getGlobalRowCopy(GO GlobalRow, std::pair<GO*,size_t> Indices, std::pair<SC*, size_t> Values, size_t &NumIndices) const {
      Teuchos::ArrayView<GO> IndicesView = Teuchos::arrayView(Indices.first, Indices.second);
      Teuchos::ArrayView<SC> ValuesView  = Teuchos::arrayView(Values.first, Values.second);

      self->getGlobalRowCopy(GlobalRow, IndicesView, ValuesView, NumIndices);
    }
}


%teuchos_rcp(Tpetra::CrsMatrix<SC,LO,GO,NO,false>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_CrsMatrix_decl.hpp"

%template(TpetraCrsMatrix) Tpetra::CrsMatrix<SC,LO,GO,NO,false>;
