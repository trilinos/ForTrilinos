/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module forteuchos

%include "copyright.i"

%include "ForTrilinosTeuchos_config.hpp"

%include <std_vector.i>
%template(VectorInt)    std::vector<int>;
%template(VectorDouble) std::vector<double>;
%template(VectorLongLong) std::vector<long long>;

// Typedefs
typedef int Teuchos_Ordinal;

%include <typemaps.i>
%fortran_view(int)
%fortran_view(long long)
%fortran_view(double)
%fortran_view(size_t)

%include "Teuchos_Types.i"
%include "Teuchos_Exceptions.i"
%include "Teuchos_RCP.i"

%include "Teuchos_Array.i"
%include "Teuchos_ArrayView.i"
%include "Teuchos_Comm.i"
%include "Teuchos_ParameterList.i"
%include "Teuchos_XML.i"
