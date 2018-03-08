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
#include "Tpetra_MultiVector.hpp"
%}

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
%ignore Tpetra::MultiVector::elementWiseMultiply;   // needs Vector
%ignore Tpetra::MultiVector::get2dCopy;             // needs ArrayView<ArrayView>
%ignore Tpetra::MultiVector::get2dView;             // needs ArrayRCP<ArrayRCP>
%ignore Tpetra::MultiVector::get2dViewNonConst;     // needs ArrayRCP<ArrayRCP>
%ignore Tpetra::MultiVector::getDualView;           // needs Kokkos::DualView
%ignore Tpetra::MultiVector::getLocalView;          // needs Kokkos::View
%ignore Tpetra::MultiVector::getVector;             // needs Tpetra::Vector
%ignore Tpetra::MultiVector::getVectorNonConst;     // needs Tpetra::Vector
%ignore Tpetra::MultiVector::modify;                // templated on device type
%ignore Tpetra::MultiVector::need_sync;             // templated on device type
%ignore Tpetra::MultiVector::scale(const Kokkos::View<const impl_scalar_type*, device_type>& alpha);
%ignore Tpetra::MultiVector::sync;                  // templated on device type
%ignore Tpetra::MultiVector::operator=;
%ignore Tpetra::MultiVector::subCopy(const Teuchos::Range1D &colRng) const; // prefer ArrayView version
%ignore Tpetra::MultiVector::subView(const Teuchos::Range1D &colRng) const; // prefer ArrayView version
%ignore Tpetra::MultiVector::subViewNonConst(const Teuchos::Range1D &colRng); // prefer ArrayView version
%ignore Tpetra::MultiVector::normWeighted;  // deprecated in Tpetra
%ignore Tpetra::deep_copy;

// =======================================================================
// Fix ±1 issues
// =======================================================================
%typemap(in)  size_t j %{$1 = *$input - 1;%}
%typemap(in)  size_t col %{$1 = *$input - 1;%}
%typemap(in)  int lclRow %{$1 = *$input - 1;%} /* int = LocalOrdinal */

// =======================================================================
// Make interface more Fortran friendly
// =======================================================================
%extend Tpetra::MultiVector<SC,LO,GO,NO> {
    MultiVector (const Teuchos::RCP< const map_type > &map, std::pair<const SC*,size_t> A, const size_t LDA, const size_t NumVectors) {
      Teuchos::ArrayView<const SC> AView = Teuchos::arrayView(A.first, A.second);
      return new Tpetra::MultiVector<SC,LO,GO,NO>(map, AView, LDA, NumVectors);
    }
    std::pair<const SC*,size_t> getData(size_t j) const {
      Teuchos::ArrayRCP<const SC> a = self->getData(j);
      return std::make_pair<const SC*,size_t>(a.get(), a.size());
    }
    std::pair<SC*,size_t> getDataNonConst(size_t j) {
      Teuchos::ArrayRCP<SC> a = self->getDataNonConst(j);
      return std::make_pair<SC*,size_t>(a.get(), a.size());
    }
    Teuchos::RCP<Tpetra::MultiVector<SC,LO,GO,NO> > subCopy(std::pair<const size_t*,size_t> cols) const {
      Teuchos::Array<size_t> colsArray(cols.second);
      for (int i = 0; i < colsArray.size(); i++)
        colsArray[i] = cols.first[i]-1;
      return self->subCopy(colsArray);
    }
    Teuchos::RCP<const Tpetra::MultiVector<SC,LO,GO,NO> > subView(std::pair<const size_t*,size_t> cols) const {
      Teuchos::Array<size_t> colsArray(cols.second);
      for (int i = 0; i < colsArray.size(); i++)
        colsArray[i] = cols.first[i]-1;
      return self->subView(colsArray);
    }
    Teuchos::RCP<Tpetra::MultiVector<SC,LO,GO,NO> > subViewNonConst(std::pair<const size_t*,size_t> cols) {
      Teuchos::Array<size_t> colsArray(cols.second);
      for (int i = 0; i < colsArray.size(); i++)
        colsArray[i] = cols.first[i]-1;
      return self->subViewNonConst(colsArray);
    }
    void dot( const Tpetra::MultiVector<SC,LO,GO,NO> &A, std::pair<SC*,size_t> dots) const {
      Teuchos::ArrayView<SC> dotsView = Teuchos::arrayView(dots.first, dots.second);
      return self->dot(A, dotsView);
    }
    void norm1(std::pair<SC*,size_t> norms) const {
      Teuchos::ArrayView<SC> normsView = Teuchos::arrayView(norms.first, norms.second);
      return self->norm1(normsView);
    }
    void norm2(std::pair<SC*,size_t> norms) const {
      Teuchos::ArrayView<SC> normsView = Teuchos::arrayView(norms.first, norms.second);
      return self->norm2(normsView);
    }
    void normInf(std::pair<SC*,size_t> norms) const {
      Teuchos::ArrayView<SC> normsView = Teuchos::arrayView(norms.first, norms.second);
      return self->normInf(normsView);
    }
    void scale(std::pair<const SC*,size_t> alpha) {
      Teuchos::ArrayView<const SC> alphaView = Teuchos::arrayView(alpha.first, alpha.second);
      self->scale(alphaView);
    }
    void meanValue(std::pair<SC*,size_t> means) const {
      Teuchos::ArrayView<SC> meansView = Teuchos::arrayView(means.first, means.second);
      self->meanValue(meansView);
    }
    void get1dCopy(std::pair<SC*,size_t> A, const size_t LDA) const {
      Teuchos::ArrayView<SC> AView = Teuchos::arrayView(A.first, A.second);
      self->get1dCopy(AView, LDA);
    }
    std::pair<const SC*,size_t> get1dView() const {
      auto a = self->get1dView();
      return std::make_pair<const SC*,size_t>(a.getRawPtr(), a.size());
    }
    std::pair<SC*,size_t> get1dViewNonConst() {
      auto a = self->get1dViewNonConst();
      return std::make_pair<SC*,size_t>(a.getRawPtr(), a.size());
    }
    void doImport (const Tpetra::MultiVector<SC,LO,GO,NO> &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      self->doImport(source, importer, CM);
    }
    void doImport (const Tpetra::MultiVector<SC,LO,GO,NO> &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      self->doImport(source, exporter, CM);
    }
    void doExport (const Tpetra::MultiVector<SC,LO,GO,NO> &source, const Tpetra::Export< LO, GO, NO > &exporter, CombineMode CM) {
      self->doExport(source, exporter, CM);
    }
    void doExport (const Tpetra::MultiVector<SC,LO,GO,NO> &source, const Tpetra::Import< LO, GO, NO > &importer, CombineMode CM) {
      self->doExport(source, importer, CM);
    }
}
%ignore Tpetra::MultiVector::MultiVector (const Teuchos::RCP< const map_type > &map, const Teuchos::ArrayView< const Scalar > &A, const size_t LDA, const size_t NumVectors);
%ignore Tpetra::MultiVector::subCopy(const Teuchos::ArrayView< const size_t > &cols) const;
%ignore Tpetra::MultiVector::subView(const Teuchos::ArrayView< const size_t > &cols) const;
%ignore Tpetra::MultiVector::subViewNonConst(const Teuchos::ArrayView< const size_t > &cols);
%ignore Tpetra::MultiVector::dot;
%ignore Tpetra::MultiVector::norm1;
%ignore Tpetra::MultiVector::norm2;
%ignore Tpetra::MultiVector::normInf;
%ignore Tpetra::MultiVector::meanValue;
%ignore Tpetra::MultiVector::scale(const Teuchos::ArrayView< const Scalar > &alpha);
%ignore Tpetra::MultiVector::getData;
%ignore Tpetra::MultiVector::getDataNonConst;
%ignore Tpetra::MultiVector::get1dCopy;
%ignore Tpetra::MultiVector::get1dView;
%ignore Tpetra::MultiVector::get1dViewNonConst;

/* Include the multivector *before* the RCP declaration so that
 * SWIG becomes aware of the default template arguments */
%include "Tpetra_MultiVector_decl.hpp"

%teuchos_rcp(Tpetra::MultiVector<SC,LO,GO,NO>)
%template(TpetraMultiVector) Tpetra::MultiVector<SC,LO,GO,NO>;
