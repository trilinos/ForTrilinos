// Dependencies
%include "Teuchos_RCP.i"
%import <Teuchos_ArrayView.i>
%import <Teuchos_Comm.i>

%{
#include "Tpetra_MultiVector.hpp"
%}

#define TPETRA_DEPRECATED

// ignore Details namespace
%ignore Details;

// POSTPONE
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> >&ArrayOfPtrs, const size_t NumVectors);
%ignore Tpetra::MultiVector::MultiVector(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>& source, const Teuchos::DataAccess copyOrView);
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, const dual_view_type& view, const Teuchos::ArrayView<const size_t>& whichVectors);
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, const dual_view_type& view, const dual_view_type& origView, const Teuchos::ArrayView<const size_t>& whichVectors);
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, const dual_view_type& view);
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, const typename dual_view_type::t_dev& d_view);
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map, const dual_view_type& view, const dual_view_type& origView);
%ignore Tpetra::MultiVector::assign;
%ignore Tpetra::MultiVector::get1dCopy;
%ignore Tpetra::MultiVector::get2dCopy;
%ignore Tpetra::MultiVector::get1dView;
%ignore Tpetra::MultiVector::get1dViewNonConst;
%ignore Tpetra::MultiVector::get2dView;
%ignore Tpetra::MultiVector::get2dViewNonConst;
%ignore Tpetra::MultiVector::offsetView;
%ignore Tpetra::MultiVector::subCopy;
%ignore Tpetra::MultiVector::subView;
%ignore Tpetra::MultiVector::subViewNonConst;
%ignore Tpetra::MultiVector::sync;
%ignore Tpetra::MultiVector::need_sync;
%ignore Tpetra::MultiVector::modify;
%ignore Tpetra::MultiVector::getLocalView;
%ignore Tpetra::MultiVector::dot;
%ignore Tpetra::MultiVector::norm1;
%ignore Tpetra::MultiVector::norm2;
%ignore Tpetra::MultiVector::normInf;
%ignore Tpetra::MultiVector::describe;
%ignore Tpetra::MultiVector::setCopyOrView;
%ignore Tpetra::MultiVector::getCopyOrView;
%ignore Tpetra::MultiVector::operator=;
%ignore Tpetra::MultiVector::getData;
%ignore Tpetra::MultiVector::getDataNonConst;
%ignore Tpetra::MultiVector::getDualView;
%ignore Tpetra::MultiVector::getVector;
%ignore Tpetra::MultiVector::getVectorNonConst;
%ignore Tpetra::MultiVector::elementWiseMultiply;
%ignore Tpetra::MultiVector::scale (const Kokkos::View<const impl_scalar_type*, device_type>& alpha);

%teuchos_rcp(Tpetra::MultiVector<SC,LO,GO,NO,false>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_MultiVector_decl.hpp"

%template(TpetraMultiVector) Tpetra::MultiVector<SC,LO,GO,NO,false>;
