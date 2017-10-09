// Dependencies
%include "Teuchos_RCP.i"
%import <Teuchos_ArrayView.i>
%import <Teuchos_Comm.i>

%{
#include "Tpetra_CrsMatrix.hpp"
%}

#define TPETRA_DEPRECATED

#define KOKKOS_INLINE_FUNCTION

%ignore Tpetra::CrsMatrix::pack_functor;

// POSTPONE
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap, \
               const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc, \
               ProfileType pftype = DynamicProfile, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::ArrayRCP<const size_t>& NumEntriesPerRowToAlloc, \
               ProfileType pftype = DynamicProfile, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap, \
               const typename local_matrix_type::row_map_type& rowPointers, \
               const typename local_graph_type::entries_type::non_const_type& columnIndices, \
               const typename local_matrix_type::values_type& values, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap, \
               const Teuchos::ArrayRCP<size_t>& rowPointers, \
               const Teuchos::ArrayRCP<LocalOrdinal>& columnIndices, \
               const Teuchos::ArrayRCP<Scalar>& values, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const crs_graph_type>& graph, \
                        const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap, \
               const local_matrix_type& lclMatrix, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
%ignore Tpetra::CrsMatrix::CrsMatrix (const Teuchos::RCP<const map_type>& rowMap, \
               const Teuchos::RCP<const map_type>& colMap, \
               const Teuchos::RCP<const map_type>& domainMap, \
               const Teuchos::RCP<const map_type>& rangeMap, \
               const local_matrix_type& lclMatrix, \
               const Teuchos::RCP<Teuchos::ParameterList>& params = Teuchos::null);
%ignore Tpetra::CrsMatrix::setAllValues (const Teuchos::ArrayRCP<size_t>& ptr, \
                  const Teuchos::ArrayRCP<LocalOrdinal>& ind, \
                  const Teuchos::ArrayRCP<Scalar>& val);
%ignore Tpetra::CrsMatrix::setAllValues(const typename local_matrix_type::row_map_type& ptr, \
                  const typename local_graph_type::entries_type::non_const_type& ind, \
                  const typename local_matrix_type::values_type& val);
%ignore Tpetra::CrsMatrix::getAllValues (Teuchos::ArrayRCP<const size_t>& rowPointers, \
                  Teuchos::ArrayRCP<const LocalOrdinal>& columnIndices, \
                  Teuchos::ArrayRCP<const Scalar>& values) const;
%ignore Tpetra::CrsMatrix::getAllValues (Teuchos::ArrayRCP<const size_t>& rowPointers, \
                  Teuchos::ArrayRCP<const LocalOrdinal>& columnIndices, \
                  Teuchos::ArrayRCP<const Scalar>& values) const;
%ignore Tpetra::CrsMatrix::operator();
%ignore Tpetra::CrsMatrix::describe;
%ignore Tpetra::CrsMatrix::copyAndPermute;
%ignore Tpetra::CrsMatrix::getLocalDiagCopy;
%ignore Tpetra::CrsMatrix::reindexColumns;
%ignore Tpetra::CrsMatrix::getNode;;
%ignore Tpetra::CrsMatrix::pack;
%ignore Tpetra::CrsMatrix::packAndPrepare;
%ignore Tpetra::CrsMatrix::packNonStatic;
%ignore Tpetra::CrsMatrix::unpackAndCombine;;
%ignore Tpetra::CrsMatrix::importAndFillComplete;
%ignore Tpetra::CrsMatrix::exportAndFillComplete;
%ignore Tpetra::CrsMatrix::expertStaticFillComplete;
%ignore Tpetra::CrsMatrix::getLocalDiagOffsets;
%ignore Tpetra::CrsMatrix::leftScale;
%ignore Tpetra::CrsMatrix::rightScale;
%ignore Tpetra::CrsMatrix::replaceDomainMapAndImporter;
%ignore Tpetra::CrsMatrix::getGraph;
%ignore Tpetra::CrsMatrix::getCrsGraph;
%ignore Tpetra::CrsMatrix::add;
%ignore Tpetra::CrsMatrix::getLocalValuesView;
%ignore Tpetra::CrsMatrix::checkSizes;
%ignore Tpetra::CrsMatrix::getLocalMatrix;
%ignore Tpetra::CrsMatrix::getLocalRowView;

%teuchos_rcp(Tpetra::CrsMatrix<SC,LO,GO,NO,false>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_CrsMatrix_decl.hpp"

%template(TpetraCrsMatrix) Tpetra::CrsMatrix<SC,LO,GO,NO,false>;
