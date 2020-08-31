/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module fortpetra

%include "fortrilinos_copyright.i"
%include "forerror/extern_forerror.i"

%import "forteuchos/forteuchos.i"

// Hide functions that use any unknown types (Kokkos, std::ostream, etc.)
%fortranonlywrapped;

// =======================================================================
// Type definition
// =======================================================================

#define HAVE_TPETRA_INST_INT_LONG_LONG
// From teuchos/kokkoscompat/src/KokkosCompat_ClassicNodeAPI_Wrapper.hpp

%ignore ForTrilinos::DefaultNodeType;
namespace ForTrilinos {
struct DefaultNodeType { 
  static const bool classic = false;
};
}

%{
#include <Kokkos_DefaultNode.hpp>
#include "fortpetra/ForTrilinos_DefaultNodeType.hpp"
%}
%inline %{
typedef double                                  SC;
typedef int                                     LO;
typedef long long                               GO;
typedef ForTrilinos::DefaultNodeType            NO;
typedef char                                    Packet;
%}

// NOTE: with the latest SWIG, these will *not* show up in downstream
// %imported modules.
%insert("fuse") {
 use, intrinsic :: iso_c_binding, only : &
   c_bool, &
   c_int, &
   c_long, &
   c_long_long, &
   c_size_t, &
   c_double, &
   scalar_type => c_double, &
   global_ordinal_type => c_long_long, &
   global_size_type => c_long, &
   size_type => c_size_t, &
   int_type => c_int, &
   mag_type => c_double, &
   norm_type => c_double
}
%insert("fdecl") {
public :: scalar_type
public :: global_ordinal_type
public :: global_size_type
public :: size_type
public :: int_type
public :: mag_type
public :: norm_type
}

// =======================================================================
// Typemaps
// =======================================================================

// Define typemap for converting Fortran indexing to C indexing
%typemap(in) int INDEX "$1 = *$input - 1;"
%typemap(out) int INDEX "$result = $1 + 1;"

%fragment("<type_traits>", "header") %{
#include <type_traits>
%}

// *Input* fortran indices
%typemap(in, fragment="<type_traits>", noblock=1) const Teuchos::ArrayView<const int>& INDEX
    ($1_basetype::value_type* tmpbegin,
     Teuchos::Array<std::remove_const<$1_basetype::value_type>::type> tmparr,
     $1_basetype tmpview)
{
  tmpbegin = static_cast<$1_basetype::value_type*>($input->data);
  tmparr.resize($input->size);
  for (int i = 0; i < tmparr.size(); i++)
    tmparr[i] = tmpbegin[i] - 1;
  tmpview = tmparr();
  $1 = &tmpview;
}

// *Input* fortran indices as an RCP
%typemap(in, fragment="<type_traits>", noblock=1) const Teuchos::ArrayRCP<const int>& INDEX
    ($1_basetype::value_type* tmpbegin,
     Teuchos::ArrayRCP<std::remove_const<$1_basetype::value_type>::type> tmparr,
     $1_basetype tmprcp)
{
  tmpbegin = static_cast<$1_basetype::value_type*>($input->data);
  tmparr.resize($input->size);
  for (int i = 0; i < tmparr.size(); i++)
    tmparr[i] = tmpbegin[i] - 1;
  tmprcp = tmparr;
  $1 = &tmprcp;
}

// Passing a *mutable* array view by const reference: increment the result
// before returning
%typemap(argout, noblock=1) const Teuchos::ArrayView<int>& INDEX
{
  for (int i = 0; i < tmpview$argnum.size(); i++)
    tmpview$argnum[i] += 1;
}

// =======================================================================
// Wrap
// =======================================================================

// All enums should be prefaced with Tpetra
%rename("Tpetra%s", %$isenumitem) "";
%rename("Tpetra%s", %$isenum)     "";
// Ignore implementation functions
%rename("$ignore", regextarget=1) "Impl$";

%include "Tpetra_ConfigDefs.i"

// Ignore indexBase, getNode (Map, CrsGraph, CrsMatrix)
%ignore *::getIndexBase;
%ignore *::getNode;

// Ignore depreceated functions (CrsGraph, CrsMatrix)
%ignore *::getGlobalNumDiags;
%ignore *::getNodeNumDiags;
%ignore *::isLowerTriangular;
%ignore *::isUpperTriangular;

// Order matters!!!
%include "Tpetra_Map.i"
%include "Tpetra_Import_Export.i"
/* %include "Tpetra_DistObject.i" */
%include "Tpetra_MultiVector.i"
/* %include "Tpetra_Vector.i" */        // needs better support for inheritance
%include "Tpetra_Operator.i"
%include "Tpetra_CrsGraph.i"
/* %include "Tpetra_RowMatrix.i" */     // needs better support for abstract classes
%include "Tpetra_CrsMatrix.i"

%include "Tpetra_InOut.i"
%include "Tpetra_MatrixMatrix.i"
