//---------------------------------*-SWIG-*----------------------------------//
/*!
 * \file   parameterlist/Teuchos_RCP.i
 * \author Seth R Johnson
 * \date   Thu Dec 08 10:40:21 2016
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
%{
#include "Teuchos_RCP.hpp"
%}

#define SWIG_SHARED_PTR_NAMESPACE Teuchos
#define shared_ptr RCP
%include <boost_shared_ptr.i>

%define %teuchos_rcp(CLASS...)
  %shared_ptr(CLASS)
%enddef

#if 0
// SNIP: Teuchos.i:RCP_DAP

%teuchos_rcp(std::basic_ostream)
%teuchos_rcp(std::ostream)
%teuchos_rcp(std::vector< int, std::allocator< int > >)
%teuchos_rcp(Teuchos::SerialDenseMatrix< int, double >)

// Enums
//%ignore Teuchos::ENull;
%import "Teuchos_ENull.hpp"
%ignore *(Teuchos::ENull);
//%import "Teuchos_RCPNode.hpp"
#endif

//---------------------------------------------------------------------------//
// end of parameterlist/Teuchos_RCP.i
//---------------------------------------------------------------------------//
