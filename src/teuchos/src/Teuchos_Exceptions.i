//---------------------------------*-SWIG-*----------------------------------//
/*!
 * \file   Teuchos_Exceptions.i
 * \author Seth R Johnson
 * \date   Tue May 02 11:26:38 2017
 * \note   Copyright (c) 2017 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
%include <std_except.i>

%{
#include "Teuchos_Exceptions.hpp"
%}

%ignore Teuchos::ExceptionBase;

namespace Teuchos {
class ExceptionBase : public std::logic_error { /* * */ };
}

%exception {
    // Make sure no unhandled exceptions exist before performing a new action
    swig::fortran_check_unhandled_exception();
    try
    {
        // Attempt the wrapped function call
        $action
    }
    catch (const std::range_error& e)
    {
        // Store a C++ exception
        SWIG_exception(SWIG_IndexError, e.what());
    }
    catch (const std::exception& e)
    {
        // Store a C++ exception
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
    catch (...)
    {
        SWIG_exception(SWIG_UnknownError, "An unknown exception occurred");
    }
}

//---------------------------------------------------------------------------//
// end of Teuchos_Exceptions.i
//---------------------------------------------------------------------------//
