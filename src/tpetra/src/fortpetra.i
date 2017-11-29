/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module fortpetra

%include "copyright.i"

%import <forteuchos.i>

%{
#include "Tpetra_ConfigDefs.hpp"

typedef double                                  SC;
typedef int                                     LO;
typedef long long                               GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
%}

typedef double                                  SC;
typedef int                                     LO;
typedef long long                               GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode NO;

%fragment("TpetraTypes", "fimports") {
 use, intrinsic :: iso_c_binding, only : &
   c_bool, &
   c_int, &
   c_long, &
   c_long_long, &
   c_size_t, &
   c_double, &
   scalar_type => c_double, &
   local_ordinal_type => c_int, &
   global_ordinal_type => c_long_long, &
   global_size_type => c_long, &
   size_type => c_size_t, &
   bool_type => c_bool, &
   int_type => c_int, &
   mag_type => c_double, &
   norm_type => c_double
}
%fragment("TpetraTypesPublic", "fpublic") {
public :: scalar_type
public :: local_ordinal_type
public :: global_ordinal_type
public :: global_size_type
public :: size_type
public :: bool_type
public :: int_type
public :: mag_type
public :: norm_type
}
// Insert fragment
%fragment("TpetraTypes");
%fragment("TpetraTypesPublic");

// Helper
namespace Kokkos {
  namespace Details {

    template<class T>
    class ArithTraits {
    public:
      typedef T val_type;
      typedef T mag_type;
    };

    template<class T>
    class InnerProductSpaceTraits {
    public:
      typedef T dot_type;
    };

  }
}
%template() Kokkos::Details::ArithTraits<SC>;
%template() Kokkos::Details::InnerProductSpaceTraits<SC>;

%ignore Teuchos::SerializationTraits;

// ignore Details namespace
%ignore Tpetra::Details;

// ignore these defines
#define TPETRA_DEPRECATED
#define KOKKOS_INLINE_FUNCTION

// ignore indexBase
%ignore getIndexBase;

// Order matters!!!
%include "Tpetra_Map.i"
%include "Tpetra_Export.i"
%include "Tpetra_Import.i"
%include "Tpetra_MultiVector.i"
/* %include "Tpetra_Vector.i" */        // needs better support for inheritance
/* %include "Tpetra_Operator.i" */      // needs to understand that Tpetra::MultiVector<SC,LO,GO,NO,false> and Tpetra::MultiVector<SC,LO,GO,NO> are the same thing
%include "Tpetra_CrsGraph.i"
/* %include "Tpetra_RowMatrix.i" */     // needs better support for abstract classes
%include "Tpetra_CrsMatrix.i"
