/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%{
#include <Tpetra_CrsMatrix.hpp>
%}

// =======================================================================
// Ignore permanently
// =======================================================================
// Ignore Details namespace
%ignore Details;
%ignore Tpetra::CrsMatrix::pack_functor;
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

// Don't warn about ignored overloads
%warnfilter(SWIGWARN_LANG_OVERLOAD_IGNORED) Tpetra::CrsMatrix::sumIntoGlobalValues;
%warnfilter(SWIGWARN_LANG_OVERLOAD_IGNORED) Tpetra::CrsMatrix::sumIntoLocalValues;
%warnfilter(SWIGWARN_LANG_OVERLOAD_IGNORED) Tpetra::CrsMatrix::CrsMatrix;

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::CrsMatrix::insertGlobalValues (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const Scalar vals[], const GlobalOrdinal inds[]); // prefer ArrayView version
%ignore Tpetra::CrsMatrix::insertLocalValues (const LocalOrdinal localRow, const LocalOrdinal numEnt, const Scalar vals[], const LocalOrdinal cols[]); // prefer ArrayView version
%ignore Tpetra::CrsMatrix::replaceGlobalValues (const GlobalOrdinal globalRow, const typename UnmanagedView< GlobalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals) const;  // needs Kokkos::UnmanagedView
%ignore Tpetra::CrsMatrix::replaceGlobalValues (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const Scalar vals[], const GlobalOrdinal cols[]) const;                                                         // prefer ArrayView variant
%ignore Tpetra::CrsMatrix::sumIntoGlobalValues (const GlobalOrdinal globalRow, const LocalOrdinal numEnt, const Scalar vals[], const GlobalOrdinal cols[], const bool atomic=useAtomicUpdatesByDefault); // prefer ArrayView variant
%ignore Tpetra::CrsMatrix::sumIntoLocalValues (const LocalOrdinal localRow, const typename UnmanagedView< LocalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals, const bool atomic=useAtomicUpdatesByDefault) const; // needs Kokkos::UnmanagedView
%ignore Tpetra::CrsMatrix::sumIntoLocalValues (const LocalOrdinal localRow, const LocalOrdinal numEnt, const Scalar vals[], const LocalOrdinal cols[], const bool atomic=useAtomicUpdatesByDefault) const; // prefer ArrayView variant
%ignore Tpetra::CrsMatrix::add;                             // needs Tpetra::RowMatrix
%ignore Tpetra::CrsMatrix::getGraph;                        // needs Tpetra::RowGraph
%ignore Tpetra::CrsMatrix::getLocalDiagCopy;                // needs Tpetra::Vector
%ignore Tpetra::CrsMatrix::getLocalMatrix;                  // needs KokkosSparse::CrsMatrix
%ignore Tpetra::CrsMatrix::getLocalValuesView;              // needs Kokkos::View
%ignore Tpetra::CrsMatrix::leftScale;                       // needs Tpetra::Vector
%ignore Tpetra::CrsMatrix::rightScale;                      // needs Tpetra::Vector
%ignore Tpetra::CrsMatrix::reorderedGaussSeidel;
%ignore Tpetra::CrsMatrix::reorderedGaussSeidelCopy;
%ignore Tpetra::importAndFillCompleteCrsMatrix;
%ignore Tpetra::exportAndFillCompleteCrsMatrix;

%ignore Tpetra::CrsMatrix::getLocalRowView;                 // ±1 issue
%ignore Tpetra::CrsMatrix::transformLocalValues;            // ±1 issue

// =======================================================================
// Fix ±1 issues
// =======================================================================
%apply int INDEX { int localRow, int LocalRow };

%apply const Teuchos::ArrayView<const int>& INDEX {
    const Teuchos::ArrayView<const LO>& indices,
    const Teuchos::ArrayView<const LO>& cols,
    Teuchos::ArrayView<const LO> cols}

%apply const Teuchos::ArrayView<int>& INDEX { const Teuchos::ArrayView<LO>& colInds, const Teuchos::ArrayView<LO>& colInds }

%apply const Teuchos::ArrayRCP<const int>& INDEX {
    const Teuchos::ArrayRCP<size_t>& rowPointers,
    const Teuchos::ArrayRCP<int>& columnIndices,
    const Teuchos::ArrayRCP<size_t>& ptr,
    const Teuchos::ArrayRCP<int>& ind
}

// Don't use atomics for value summation
%typemap(in,numinputs=0) bool atomic {
  $1 = false;
}

// =======================================================================
// Make interface more Fortran friendly
// =======================================================================
%extend Tpetra::CrsMatrix<SC,LO,GO,NO> {
    void getAllValues(Teuchos::ArrayView<size_t> rowPointers, Teuchos::ArrayView<LO> columnIndices, Teuchos::ArrayView<SC> values) const {
        auto rowptr = $self->getLocalRowPtrsHost();
        auto colidx = $self->getLocalIndicesHost();
        auto val = $self->getLocalValuesHost(Tpetra::Access::ReadOnly);

        TEUCHOS_TEST_FOR_EXCEPTION(rowptr.size() != rowPointers.size(),   std::runtime_error, "Wrong rowPointers size");
        TEUCHOS_TEST_FOR_EXCEPTION(colidx.size() != columnIndices.size(), std::runtime_error, "Wrong columnIndices size");
        TEUCHOS_TEST_FOR_EXCEPTION(val.size()    != values.size(),        std::runtime_error, "Wrong values size");
        auto n = rowPointers.size();
        for (int i = 0; i < n; i++)
            rowPointers[i] = rowptr[i]+1;
        auto nnz = columnIndices.size();
        for (int i = 0; i < nnz; i++) {
            columnIndices[i] = colidx[i]+1;
            values       [i] = val[i];
        }
    }
    // NOTE: This is semantically different function from Tpetra. Here, we *require* that user already allocated the arrays to store the data
    void
    getGlobalRowCopy (GO GlobalRow,
                      const Teuchos::ArrayView<GO>& Indices,
                      const Teuchos::ArrayView<SC>& Values,
                      size_t& NumEntries) const {
        Tpetra::CrsMatrix<SC,LO,GO,NO>::nonconst_global_inds_host_view_type temp_i(Indices.data(), Indices.size());
        Tpetra::CrsMatrix<SC,LO,GO,NO>::nonconst_values_host_view_type temp_v(Values.data(), Values.size());
        return $self->getGlobalRowCopy(GlobalRow, temp_i, temp_v, NumEntries);
    }

    void
    getLocalRowCopy (LO localRow,
                     const Teuchos::ArrayView<LO>& colInds,
                     const Teuchos::ArrayView<SC>& vals,
                     size_t& numEntries) const {
        Tpetra::CrsMatrix<SC,LO,GO,NO>::nonconst_local_inds_host_view_type temp_i(colInds.data(), colInds.size());
        Tpetra::CrsMatrix<SC,LO,GO,NO>::nonconst_values_host_view_type temp_v(vals.data(), vals.size());
        return $self->getLocalRowCopy(localRow, temp_i, temp_v, numEntries);
    }
}

// Add doImport and doExport
%tpetra_extend_with_import_export(Tpetra::CrsGraph<LO,GO,NO>)

%ignore Tpetra::CrsMatrix::getAllValues(Teuchos::ArrayRCP< const size_t > &rowPointers, Teuchos::ArrayRCP< const LocalOrdinal > &columnIndices, Teuchos::ArrayRCP< const Scalar > &values) const;
%ignore Tpetra::CrsMatrix::getLocalRowViewRaw(const LocalOrdinal lclRow, LocalOrdinal &numEnt, const LocalOrdinal *&lclColInds, const Scalar *&vals) const;
%ignore Tpetra::CrsMatrix::getAllValues (Teuchos::ArrayRCP<const size_t>& rowPointers, Teuchos::ArrayRCP<const LocalOrdinal>& columnIndices, Teuchos::ArrayRCP<const Scalar>& values) const;
%ignore Tpetra::CrsMatrix::sumIntoLocalValues (const LocalOrdinal localRow, const typename UnmanagedView< LocalIndicesViewType >::type &inputInds, const typename UnmanagedView< ImplScalarViewType >::type &inputVals, const bool atomic=useAtomicUpdatesByDefault) const;
%ignore Tpetra::CrsMatrix::sumIntoLocalValues (const LocalOrdinal localRow, const LocalOrdinal numEnt, const Scalar vals[], const LocalOrdinal cols[], const bool atomic=useAtomicUpdatesByDefault);


// Define type aliases for Tpetra matrix base class
namespace Tpetra {
  template <class Scalar,
            class LocalOrdinal,
            class GlobalOrdinal,
            class Node>
  class RowMatrix
{
public:
    typedef Scalar        scalar_type;
    typedef LocalOrdinal  local_ordinal_type;
    typedef GlobalOrdinal global_ordinal_type;
    typedef Node          node_type;

    using impl_scalar_type = typename Kokkos::ArithTraits<Scalar>::val_type;
    using mag_type = typename Kokkos::ArithTraits<Scalar>::mag_type;
    typedef typename
        Kokkos::View<impl_scalar_type*, typename Node::device_type>::const_type
        values_device_view_type;
    typedef typename values_device_view_type::HostMirror::const_type
        values_host_view_type;
    typedef typename values_device_view_type::HostMirror
        nonconst_values_host_view_type;

    typedef typename
        Kokkos::View<LocalOrdinal *, typename Node::device_type>::const_type
        local_inds_device_view_type;
    typedef typename local_inds_device_view_type::HostMirror::const_type
        local_inds_host_view_type;
    typedef typename local_inds_device_view_type::HostMirror
        nonconst_local_inds_host_view_type;

    typedef typename
        Kokkos::View<GlobalOrdinal *, typename Node::device_type>::const_type
        global_inds_device_view_type;
    typedef typename global_inds_device_view_type::HostMirror::const_type
        global_inds_host_view_type;
    typedef typename global_inds_device_view_type::HostMirror
        nonconst_global_inds_host_view_type;


    typedef typename
        Kokkos::View<const size_t*, typename Node::device_type>::const_type
        row_ptrs_device_view_type;
    typedef typename row_ptrs_device_view_type::HostMirror::const_type
        row_ptrs_host_view_type;
};
}

%include "Tpetra_CrsMatrix_decl.hpp"

%teuchos_rcp(Tpetra::CrsMatrix<SC,LO,GO,NO>)
%template() Tpetra::RowMatrix<SC,LO,GO,NO>;
%template(TpetraCrsMatrix) Tpetra::CrsMatrix<SC,LO,GO,NO>;

// Operator to Matrix conversion
%{
#include <Tpetra_CrsMatrix.hpp>
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
