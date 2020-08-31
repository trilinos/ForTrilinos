/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module forteuchos

%include "fortrilinos_copyright.i"
%include "forerror/extern_forerror.i"

// Convert all std::string references/values to and from Fortran strings
%include <std_string.i>

// =======================================================================
// Global treatments
// =======================================================================

// By default, wrap all constants as Fortran literals
%fortranconst;

%{
#include "ForTrilinos_config.h"
%}
%include "ForTrilinos_config.h"

// Typedefs
typedef int Teuchos_Ordinal;

// Treat plain 'int' (and Teuchos_Ordinal) as native integers.
// Assume that large integers are defined as size_t or global_ordinal
%typemap(ftype, in="integer, intent(in)") int "integer"
%typemap(fin) int "$1 = int($input, C_INT)"
%typemap(fout) int "$result = int($1)"

// Ignore methods that require non-implemented types or operators
%ignore *::describe;
%ignore *::print(std::ostream& os);
%ignore *::operator=;
%ignore *::operator<<;
%ignore *::operator();

// =======================================================================
// Wrap
// =======================================================================

// Prefix all enums with Teuchos
%rename("Teuchos%s",%$isenumitem) "";
%rename("Teuchos%s",%$isenum) "";

%include "Teuchos_Traits.i"
%include "Teuchos_Types.i"
%include "Teuchos_Exceptions.i"
%include "Teuchos_RCP.i"

%include "Teuchos_ArrayView.i"
%include "Teuchos_ArrayRCP.i"
%include "Teuchos_Array.i"
%include "Teuchos_Comm.i"
%include "Teuchos_ParameterList.i"
%include "Teuchos_XML.i"
