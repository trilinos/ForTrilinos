//---------------------------------*-SWIG-*----------------------------------//
/*!
 * \file   simple/eigen_handle.i
 * \author Seth R Johnson
 * \date   Mon Feb 27 10:43:39 2017
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
%{
#include "eigen_handle.hpp"
%}

%ignore setup_matrix(Teuchos::RCP<Matrix>);
%ignore setup_matrix_rhs(Teuchos::RCP<Matrix>);

%include "Teuchos_Comm.i"
%include "eigen_handle.hpp"

//---------------------------------------------------------------------------//
// end of simple/eigen_handle.i
//---------------------------------------------------------------------------//
