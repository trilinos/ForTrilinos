/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module fortpetra

%include <copyright.i>
%include <extern_forerror.i>

%import <forteuchos.i>

%{
#include "Tpetra_ConfigDefs.hpp"

typedef double                                  SC;
typedef int                                     LO;
typedef long long                               GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
typedef char                                    Packet;
%}

#define HAVE_TPETRA_INST_INT_LONG_LONG
typedef double                                  SC;
typedef int                                     LO;
typedef long long                               GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
typedef char                                    Packet;


// FIXME: Restore previous bool behaviour
FORT_FUND_TYPEMAP(bool, "logical(C_BOOL)")

%fragment("TpetraTypes", "fmodule") {
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

// enum workaround
#define RENAME_ENUM(X) %rename(Tpetra##X) X;
RENAME_ENUM(LocalGlobal)
RENAME_ENUM(LocallyReplicated)
RENAME_ENUM(GloballyDistributed)
RENAME_ENUM(LookupStatus)
RENAME_ENUM(AllIDsPresent)
RENAME_ENUM(IDNotPresent)
RENAME_ENUM(ProfileType)
RENAME_ENUM(StaticProfile)
RENAME_ENUM(DynamicProfile)
RENAME_ENUM(OptimizeOption)
RENAME_ENUM(DoOptimizeStorage)
RENAME_ENUM(DoNotOptimizeStorage)
RENAME_ENUM(ESweepDirection)
RENAME_ENUM(Forward)
RENAME_ENUM(Backward)
RENAME_ENUM(Symmetric)
RENAME_ENUM(CombineMode)
RENAME_ENUM(ADD)
RENAME_ENUM(INSERT)
RENAME_ENUM(REPLACE)
RENAME_ENUM(ABSMAX)
RENAME_ENUM(ZERO)
RENAME_ENUM(ELocalGlobal)
RENAME_ENUM(LocalIndices)
RENAME_ENUM(GlobalIndices)
#undef RENAME_ENUM

// ignore Details namespace
%ignore Tpetra::Details;

// ignore these defines
#define TPETRA_DEPRECATED
#define KOKKOS_INLINE_FUNCTION

// ignore indexBase
%ignore getIndexBase;

// Some enums
%ignore EPrivateComputeViewConstructor;
%ignore EPrivateHostViewConstructor;
%include "Tpetra_ConfigDefs.hpp"
%include "Tpetra_CombineMode.hpp"

// Order matters!!!
%include "Tpetra_Map.i"
%include "Tpetra_Export.i"
%include "Tpetra_Import.i"
/* %include "Tpetra_DistObject.i" */
%include "Tpetra_MultiVector.i"
/* %include "Tpetra_Vector.i" */        // needs better support for inheritance
/* %include "Tpetra_Operator.i" */      // needs to understand that Tpetra::MultiVector<SC,LO,GO,NO,false> and Tpetra::MultiVector<SC,LO,GO,NO> are the same thing
%include "Tpetra_CrsGraph.i"
/* %include "Tpetra_RowMatrix.i" */     // needs better support for abstract classes
%include "Tpetra_CrsMatrix.i"

%include "Tpetra_InOut.i"
%include "Tpetra_MatrixMatrix.i"
