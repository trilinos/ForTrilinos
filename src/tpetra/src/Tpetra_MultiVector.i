// Dependencies
%include "Teuchos_RCP.i"
%import <Teuchos_ArrayView.i>
%import <Teuchos_Comm.i>

%{
#include "Tpetra_MultiVector.hpp"
%}

// =======================================================================
// Ignore permanently
// =======================================================================
%ignore Teuchos::TypeNameTraits;

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, \
        const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> >&ArrayOfPtrs, \
        const size_t NumVectors);                               // needs ArrayView of ArrayView (or alternative)
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, \
        const dual_view_type& view, const Teuchos::ArrayView<const size_t>& whichVectors);  // needs Kokkos::DualView; needs Teuchos::ArrayView<size_t>
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, \
        const dual_view_type& view, \
        const dual_view_type& origView, \
        const Teuchos::ArrayView<const size_t>& whichVectors);  // needs Kokkos::DualView; needs Teuchos::ArrayView<size_t>
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, \
        const dual_view_type& view);                            // needs Kokkos::DualView
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, \
        const typename dual_view_type::t_dev& d_view);          // needs Kokkos::DualView
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, \
        const dual_view_type& view, const dual_view_type& origView);    // needs Kokkos::DualView
%ignore Tpetra::MultiVector::assign;
%ignore Tpetra::MultiVector::describe;              // needs Teuchos::FancyOStream
/* %ignore Tpetra::MultiVector::dot(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>& A, const Teuchos::ArrayView<dot_type>& dots) const; */
%ignore Tpetra::MultiVector::dot(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, std::vector< T > &dots) const;
%ignore Tpetra::MultiVector::dot(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Kokkos::View< dot_type *, device_type > &dots) const;
%ignore Tpetra::MultiVector::elementWiseMultiply;   // needs Vector
%ignore Tpetra::MultiVector::get1dView;             // needs Teuchos::ArrayRCP
%ignore Tpetra::MultiVector::get1dViewNonConst;     // needs Teuchos::ArrayRCP
%ignore Tpetra::MultiVector::get2dCopy;             // needs ArrayView<ArrayView>
%ignore Tpetra::MultiVector::get2dView;             // needs ArrayRCP<ArrayRCP>
%ignore Tpetra::MultiVector::get2dViewNonConst;     // needs ArrayRCP<ArrayRCP>
%ignore Tpetra::MultiVector::getData;               // needs Teuchos::ArrayRCP
%ignore Tpetra::MultiVector::getDataNonConst;       // needs Teuchos::ArrayRCP
%ignore Tpetra::MultiVector::getDualView;           // needs Kokkos::DualView
%ignore Tpetra::MultiVector::getLocalView;          // needs Kokkos::View
%ignore Tpetra::MultiVector::getVector;             // needs Tpetra::Vector
%ignore Tpetra::MultiVector::getVectorNonConst;     // needs Tpetra::Vector
%ignore Tpetra::MultiVector::modify;                // templated on device type
%ignore Tpetra::MultiVector::need_sync;             // templated on device type
%ignore Tpetra::MultiVector::sync;                  // templated on device type
%ignore Tpetra::MultiVector::norm1(const Kokkos::View<mag_type*, device_type>& norms) const;
%ignore Tpetra::MultiVector::norm1(const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;
%ignore Tpetra::MultiVector::norm1(const Kokkos::View<T*, device_type>& norms) const;
%ignore Tpetra::MultiVector::norm1(const Teuchos::ArrayView<T>& norms) const;
%ignore Tpetra::MultiVector::norm2(const Kokkos::View<mag_type*, device_type>& norms) const;
%ignore Tpetra::MultiVector::norm2(const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;
%ignore Tpetra::MultiVector::norm2(const Kokkos::View<T*, device_type>& norms) const;
%ignore Tpetra::MultiVector::norm2(const Teuchos::ArrayView<T>& norms) const;
%ignore Tpetra::MultiVector::normInf(const Kokkos::View<mag_type*, device_type>& norms) const;
%ignore Tpetra::MultiVector::normInf(const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;
%ignore Tpetra::MultiVector::normInf(const Kokkos::View<T*, device_type>& norms) const;
%ignore Tpetra::MultiVector::normInf(const Teuchos::ArrayView<T>& norms) const;
%ignore Tpetra::MultiVector::offsetView;
%ignore Tpetra::MultiVector::operator=;
%ignore Tpetra::MultiVector::scale(const Kokkos::View<const impl_scalar_type*, device_type>& alpha);
%ignore Tpetra::MultiVector::subCopy;               // ±1 issue; needs Teuchos::Range1D
%ignore Tpetra::MultiVector::subView;               // ±1 issue; needs Teuchos::Range1D
%ignore Tpetra::MultiVector::subViewNonConst;       // ±1 issue; needs Teuchos::Range1D

// =======================================================================
// Fix ±1 issues
// =======================================================================
%typemap(in)  const size_t j %{$1 = *$input - 1;%}

%teuchos_rcp(Tpetra::MultiVector<SC,LO,GO,NO,false>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_MultiVector_decl.hpp"

%template(TpetraMultiVector) Tpetra::MultiVector<SC,LO,GO,NO,false>;
