/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include <Teuchos_Exceptions.hpp>
%}

%ignore Teuchos::ExceptionBase;

namespace Teuchos {
class ExceptionBase : public std::logic_error { /* * */ };
}

%exception {
    // Make sure no unhandled exceptions exist before performing a new action
    SWIG_check_unhandled_exception();
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
