/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include "Teuchos_RCP.hpp"
%}

#define SWIG_SHARED_PTR_NAMESPACE Teuchos
#define shared_ptr RCP
#define SWIG_SHARED_PTR_NOT_NULL(f) (!Teuchos::is_null(f))

// Override null deleters defined in shared_ptr.i, used for storing
// references/raw pointers in shared_ptr.
%fragment("SWIG_null_deleter", "header") {
%#define SWIG_NO_NULL_DELETER_0 , Teuchos::RCP_WEAK_NO_DEALLOC
%#define SWIG_NO_NULL_DELETER_1
%#define SWIG_NO_NULL_DELETER_SWIG_POINTER_NEW
%#define SWIG_NO_NULL_DELETER_SWIG_POINTER_OWN
}

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
