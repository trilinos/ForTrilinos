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

#define HAVE_TPETRA_INST_INT_LONG_LONG
// From teuchos/kokkoscompat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp
%ignore KokkosSerialWrapperNode;
namespace Kokkos { namespace Compat {
struct KokkosSerialWrapperNode {
  static const bool classic = false;
};
} }

%{
#include "Kokkos_DefaultNode.hpp"
%}
%inline %{
typedef double                                  SC;
typedef int                                     LO;
typedef long long                               GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
typedef char                                    Packet;
%}

// NOTE: with the latest SWIG, these will *not* show up in downstream
// %imported modules even if %insert is replaced with %fragment.
%insert("fuse") {
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
   int_type => c_int, &
   mag_type => c_double, &
   norm_type => c_double
}
%insert("fdecl") {
public :: scalar_type
public :: local_ordinal_type
public :: global_ordinal_type
public :: global_size_type
public :: size_type
public :: int_type
public :: mag_type
public :: norm_type
}

// All enums should be prefaced with Tpetra
%rename("Tpetra%s", %$isenumitem) "";
%rename("Tpetra%s", %$isenum)     "";

%include "Tpetra_ConfigDefs.i"

// ignore indexBase (Map, CrsGraph, CrsMatrix)
%ignore getIndexBase;

// Order matters!!!
%include "Tpetra_Map.i"
%include "Tpetra_Export.i"
%include "Tpetra_Import.i"
/* %include "Tpetra_DistObject.i" */
%include "Tpetra_MultiVector.i"
/* %include "Tpetra_Vector.i" */        // needs better support for inheritance
%include "Tpetra_Operator.i"
%include "Tpetra_CrsGraph.i"
/* %include "Tpetra_RowMatrix.i" */     // needs better support for abstract classes
%include "Tpetra_CrsMatrix.i"

%include "Tpetra_InOut.i"
%include "Tpetra_MatrixMatrix.i"
