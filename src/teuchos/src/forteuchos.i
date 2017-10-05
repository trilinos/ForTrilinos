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

%{
#include "Teuchos_BLAS_types.hpp"
%}
%include "Teuchos_BLAS_types.hpp"

%include "Teuchos_Array.i"
%include "Teuchos_ArrayView.i"
%include "Teuchos_Comm.i"
%include "Teuchos_Exceptions.i"
%include "Teuchos_RCP.i"
%include "Teuchos_ParameterList.i"
%include "Teuchos_XML.i"

//---------------------------------------------------------------------------//
// end of parameterlist/teuchos.i
//---------------------------------------------------------------------------//
