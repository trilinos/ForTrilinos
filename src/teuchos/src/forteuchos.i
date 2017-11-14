//---------------------------------*-SWIG-*----------------------------------//
/*!
 * \file   parameterlist/forteuchos.i
 * \author Seth R Johnson
 * \date   Tue Dec 06 17:54:39 2016
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//

%module forteuchos

%include "ForTrilinosTeuchos_config.hpp"

%include <std_vector.i>
%template(VectorInt)    std::vector<int>;
%template(VectorDouble) std::vector<double>;

// Typedefs
typedef int Teuchos_Ordinal;

%include <typemaps.i>
%fortran_view(int)
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

//---------------------------------------------------------------------------//
// end of parameterlist/teuchos.i
//---------------------------------------------------------------------------//
