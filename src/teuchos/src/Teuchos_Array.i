/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include "Teuchos_Array.hpp"
%}

// Include typemap for const array references
%include <std_container.i>

namespace Teuchos
{
template<typename T>
class Array
{
    public:
    typedef T value_type;
    typedef std::size_t size_type;
};
}

%define TEUCHOS_ARRAY(TYPE...)

// Add native wrapping typemaps to convert to/from Teuchos array
%std_native_container(Teuchos::Array<TYPE>)
// Instantiate the typemaps without generating wrappers
%template() Teuchos::Array<TYPE>;

%enddef

TEUCHOS_ARRAY(int);
TEUCHOS_ARRAY(double);
TEUCHOS_ARRAY(long long);
