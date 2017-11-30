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
%ignore Tpetra::MultiVector::operator=;
%ignore Tpetra::MultiVector::scale(const Kokkos::View<const impl_scalar_type*, device_type>& alpha);    // needs Kokkos::View
%ignore Tpetra::MultiVector::subCopy(const Teuchos::Range1D &colRng) const; // prefer ArrayView version
%ignore Tpetra::MultiVector::subView(const Teuchos::Range1D &colRng) const; // prefer ArrayView version
%ignore Tpetra::MultiVector::subViewNonConst(const Teuchos::Range1D &colRng); // prefer ArrayView version
%ignore Tpetra::MultiVector::normWeighted;  // deprecated in Tpetra

// =======================================================================
// Fix ±1 issues
// =======================================================================
%typemap(in)  const size_t j %{$1 = *$input - 1;%}

// =======================================================================
// Make interface more Fortran friendly
// =======================================================================
%ignore Tpetra::MultiVector::MultiVector (const Teuchos::RCP< const map_type > &map, const Teuchos::ArrayView< const Scalar > &A, const size_t LDA, const size_t NumVectors);
%ignore Tpetra::MultiVector::subCopy(const Teuchos::ArrayView< const size_t > &cols) const;
%ignore Tpetra::MultiVector::subView(const Teuchos::ArrayView< const size_t > &cols) const;
%ignore Tpetra::MultiVector::subViewNonConst(const Teuchos::ArrayView< const size_t > &cols);
%ignore Tpetra::MultiVector::dot(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node, classic>& A, const Teuchos::ArrayView<dot_type>& dots) const;
%ignore Tpetra::MultiVector::dot(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, std::vector< T > &dots) const;
%ignore Tpetra::MultiVector::dot(const MultiVector< Scalar, LocalOrdinal, GlobalOrdinal, Node, classic > &A, const Kokkos::View< dot_type *, device_type > &dots) const;
%ignore Tpetra::MultiVector::norm1(const Kokkos::View<mag_type*, device_type>& norms) const;
%ignore Tpetra::MultiVector::norm1(const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;
%ignore Tpetra::MultiVector::norm1(const Kokkos::View<T*, device_type>& norms) const;
%ignore Tpetra::MultiVector::norm1(const Teuchos::ArrayView<mag_type>& norms) const;
%ignore Tpetra::MultiVector::norm1(const Teuchos::ArrayView<T>& norms) const;
%ignore Tpetra::MultiVector::norm2(const Kokkos::View<mag_type*, device_type>& norms) const;
%ignore Tpetra::MultiVector::norm2(const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;
%ignore Tpetra::MultiVector::norm2(const Kokkos::View<T*, device_type>& norms) const;
%ignore Tpetra::MultiVector::norm2(const Teuchos::ArrayView<mag_type>& norms) const;
%ignore Tpetra::MultiVector::norm2(const Teuchos::ArrayView<T>& norms) const;
%ignore Tpetra::MultiVector::normInf(const Kokkos::View<mag_type*, device_type>& norms) const;
%ignore Tpetra::MultiVector::normInf(const Kokkos::View<mag_type*, Kokkos::HostSpace>& norms) const;
%ignore Tpetra::MultiVector::normInf(const Kokkos::View<T*, device_type>& norms) const;
%ignore Tpetra::MultiVector::normInf(const Teuchos::ArrayView<mag_type>& norms) const;
%ignore Tpetra::MultiVector::normInf(const Teuchos::ArrayView<T>& norms) const;
%ignore Tpetra::MultiVector::meanValue(const Teuchos::ArrayView< impl_scalar_type > &means) const;
%ignore Tpetra::MultiVector::scale(const Teuchos::ArrayView< const Scalar > &alpha);
%extend Tpetra::MultiVector<SC,LO,GO,NO,false> {
    MultiVector (const Teuchos::RCP< const map_type > &map, std::pair<const SC*,size_t> A, const size_t LDA, const size_t NumVectors) {
      Teuchos::ArrayView<const SC> AView = Teuchos::arrayView(A.first, A.second);
      return new Tpetra::MultiVector<SC,LO,GO,NO,false>(map, AView, LDA, NumVectors);
    }
    Teuchos::RCP<Tpetra::MultiVector<SC,LO,GO,NO,false> > subCopy(std::pair<const size_t*,size_t> cols) const {
      Teuchos::Array<size_t> colsArray(cols.second);
      for (int i = 0; i < colsArray.size(); i++)
        colsArray[i] = cols.first[i]-1;
      return self->subCopy(colsArray);
    }
    Teuchos::RCP<const Tpetra::MultiVector<SC,LO,GO,NO,false> > subView(std::pair<const size_t*,size_t> cols) const {
      Teuchos::Array<size_t> colsArray(cols.second);
      for (int i = 0; i < colsArray.size(); i++)
        colsArray[i] = cols.first[i]-1;
      return self->subView(colsArray);
    }
    Teuchos::RCP<Tpetra::MultiVector<SC,LO,GO,NO,false> > subViewNonConst(std::pair<const size_t*,size_t> cols) {
      Teuchos::Array<size_t> colsArray(cols.second);
      for (int i = 0; i < colsArray.size(); i++)
        colsArray[i] = cols.first[i]-1;
      return self->subViewNonConst(colsArray);
    }
    void dot( const Tpetra::MultiVector<SC,LO,GO,NO,false> &A, std::pair<SC*,size_t> dots) const {
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
}
%typemap(in)  const size_t j %{$1 = *$input - 1;%}
%typemap(in)  const size_t col %{$1 = *$input - 1;%}
%typemap(in)  const int lclRow %{$1 = *$input - 1;%} /* int = LocalOrdinal */


%teuchos_rcp(Tpetra::MultiVector<SC,LO,GO,NO,false>)

#define HAVE_TPETRA_INST_INT_INT
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_MultiVector_decl.hpp"

%template(TpetraMultiVector) Tpetra::MultiVector<SC,LO,GO,NO,false>;
