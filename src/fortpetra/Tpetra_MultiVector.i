/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

%{
#include <Tpetra_MultiVector.hpp>
%}

// Treat array RCP return values as array views
// (they will not reference count, though, of course)
%apply Teuchos::ArrayView<double> { Teuchos::ArrayRCP<double> };
%apply Teuchos::ArrayView<const double> { Teuchos::ArrayRCP<const double> };

// Function signatures for local quantities are incorrectly declared as size_t
%apply LO { size_t getLocalLength,
            size_t getOrigNumLocalRows,
            size_t getOrigNumLocalCols };

// =======================================================================
// Postpone temporarily
// =======================================================================
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map,
        const Teuchos::ArrayView<const Teuchos::ArrayView<const Scalar> >&ArrayOfPtrs,
        const size_t NumVectors);                               // needs ArrayView of ArrayView (or alternative)
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map,
        const dual_view_type& view, const Teuchos::ArrayView<const size_t>& whichVectors);  // needs Kokkos::DualView; needs Teuchos::ArrayView<size_t>
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map,
        const dual_view_type& view,
        const dual_view_type& origView,
        const Teuchos::ArrayView<const size_t>& whichVectors);  // needs Kokkos::DualView; needs Teuchos::ArrayView<size_t>
%ignore Tpetra::MultiVector::MultiVector(const MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
        const map_type& subMap,
        const size_t offset = 0);                               // prefer RCP version
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map,
        const dual_view_type& view);                            // needs Kokkos::DualView
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map,
        const typename dual_view_type::t_dev& d_view);          // needs Kokkos::DualView
%ignore Tpetra::MultiVector::MultiVector(const Teuchos::RCP<const map_type>& map,
        const dual_view_type& view, const dual_view_type& origView);    // needs Kokkos::DualView
%ignore Tpetra::MultiVector::MultiVector(MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&&);     // move constructor
%ignore Tpetra::MultiVector::operator=(MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>&&);       // move assignment
%ignore Tpetra::MultiVector::assign;
%ignore Tpetra::MultiVector::describe;              // needs Teuchos::FancyOStream
%ignore Tpetra::MultiVector::elementWiseMultiply;   // needs Vector
%ignore Tpetra::MultiVector::get2dCopy;             // needs ArrayView<ArrayView>
%ignore Tpetra::MultiVector::get2dView;             // needs ArrayRCP<ArrayRCP>
%ignore Tpetra::MultiVector::get2dViewNonConst;     // needs ArrayRCP<ArrayRCP>
%ignore Tpetra::MultiVector::getDualView;           // needs Kokkos::DualView
%ignore Tpetra::MultiVector::getLocalView;          // needs Kokkos::View
%ignore Tpetra::MultiVector::getLocalViewHost;      // needs Kokkos::View
%ignore Tpetra::MultiVector::getLocalViewDevice;    // needs Kokkos::View
%ignore Tpetra::MultiVector::getVector;             // needs Tpetra::Vector
%ignore Tpetra::MultiVector::getVectorNonConst;     // needs Tpetra::Vector
%ignore Tpetra::MultiVector::modify;                // templated on device type
%ignore Tpetra::MultiVector::need_sync;             // templated on device type
%ignore Tpetra::MultiVector::sync_host;             // not needed
%ignore Tpetra::MultiVector::sync_device;           // not needed
%ignore Tpetra::MultiVector::need_sync_host;        // not needed
%ignore Tpetra::MultiVector::need_sync_device;      // not needed
%ignore Tpetra::MultiVector::modify_device;         // not needed
%ignore Tpetra::MultiVector::modify_host;           // not needed
%ignore Tpetra::MultiVector::scale(const Kokkos::View<const impl_scalar_type*, device_type>& alpha);
%ignore Tpetra::MultiVector::sync;                  // templated on device type
%ignore Tpetra::MultiVector::operator=;
%ignore Tpetra::MultiVector::subCopy(const Teuchos::Range1D &colRng) const; // prefer ArrayView version
%ignore Tpetra::MultiVector::subView(const Teuchos::Range1D &colRng) const; // prefer ArrayView version
%ignore Tpetra::MultiVector::subViewNonConst(const Teuchos::Range1D &colRng); // prefer ArrayView version
%ignore Tpetra::MultiVector::normWeighted;  // deprecated in Tpetra
%ignore Tpetra::deep_copy;
%ignore Tpetra::getMultiVectorWhichVectors;

// Add doImport and doExport
%tpetra_extend_with_import_export(Tpetra::MultiVector<SC,LO,GO,NO>)

// Fix ±1 issues
%apply int INDEX { size_t j, size_t col, int lclRow }
%apply const Teuchos::ArrayView<const int>& INDEX { const Teuchos::ArrayView<const size_t>& cols }

/* Include the multivector *before* the RCP declaration so that
 * SWIG becomes aware of the default template arguments */
%include "Tpetra_MultiVector_decl.hpp"

%teuchos_rcp(Tpetra::MultiVector<SC,LO,GO,NO>)
%template(TpetraMultiVector) Tpetra::MultiVector<SC,LO,GO,NO>;

%clear Teuchos::ArrayRCP<double>;
%clear Teuchos::ArrayRCP<const double>;
