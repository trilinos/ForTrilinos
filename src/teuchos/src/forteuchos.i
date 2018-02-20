/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module forteuchos

%include <copyright.i>
%include <extern_forerror.i>

%include "ForTrilinosTeuchos_config.hpp"

// Typedefs
typedef int Teuchos_Ordinal;

// pair<int> typemaps: (TODO: these should be deprecated)
%include <typemaps.i>
%fortran_view(int)
%fortran_view(long long)
%fortran_view(double)
%fortran_view(size_t)

// FIXME: Restore previous bool behaviour
FORT_FUND_TYPEMAP(bool, "logical(C_BOOL)")

// Convert all std::string references/values to and from Fortran strings
%include <std_string.i>

// Declare and ignore some traits classes
namespace Teuchos {
  template<typename O, typename T> class SerializationTraits;
  template<typename T> class TypeNameTraits;
}
%ignore Teuchos::SerializationTraits;
%ignore Teuchos::TypeNameTraits;

// Prefix all enums with Teuchos
%rename("Teuchos%s",%$isenumitem) "";
%rename("Teuchos%s",%$isenum) "";

%include "Teuchos_Types.i"
%include "Teuchos_Exceptions.i"
%include "Teuchos_RCP.i"

%include "Teuchos_Array.i"
%include "Teuchos_ArrayView.i"
%include "Teuchos_Comm.i"
%include "Teuchos_ParameterList.i"
%include "Teuchos_XML.i"
