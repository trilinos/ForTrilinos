/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */

// Dependencies
%include "Teuchos_RCP.i"

%{
#include "Teuchos_Array.hpp"
%}

// Typedefs
typedef int Teuchos_Ordinal;

namespace Teuchos {
  template<typename T>
  class TypeNameTraits;
}

// Make RCPs
%teuchos_rcp(Teuchos::Array<int>)
%teuchos_rcp(Teuchos::Array<double>)

// Ignore misc classes
%ignore Teuchos::TypeNameTraits;
%ignore Teuchos::InvalidArrayStringRepresentation;

// Ignore functions
%ignore operator>>;
%ignore extractDataFromISS;
%ignore getArrayTypeNameTraitsFormat;
%ignore Array(InputIterator first, InputIterator last);
%ignore Teuchos::Array::assign;
%ignore Teuchos::Array::begin;
%ignore Teuchos::Array::end;
%ignore Teuchos::Array::rbegin;
%ignore Teuchos::Array::rend;
%ignore Teuchos::Array::insert;
%ignore Teuchos::Array::erase;
%ignore Teuchos::Array::getRawPtr;
%ignore Teuchos::Array::hasBoundsChecking;
%ignore Teuchos::Array::toVector;

// Postpone
%ignore Array(const ArrayView<const T>& a);
%ignore Teuchos::Array::operator=(const Array<T>& a);
%ignore Teuchos::Array::operator=( const std::vector<T> &v );
%ignore Teuchos::Array::operator ArrayView<T>();
%ignore Teuchos::Array::operator ArrayView<const T>() const;
%ignore Teuchos::Array::view;
%ignore Teuchos::Array::operator[];
%ignore Teuchos::Array::operator();
%ignore Teuchos::Array::at;
%ignore Teuchos::Array::back;
%ignore Teuchos::Array::front;
%ignore Teuchos::Array::append;
%ignore Teuchos::Array::remove;
%ignore Teuchos::Array::toString;

%define ADD_VIEW(TYPE)
// Extend ArrayView
%extend Teuchos::Array<TYPE> {
Array(std::pair<TYPE*, size_t> view)
{
  Teuchos::Array<TYPE>* arr = new Teuchos::Array<TYPE>(view.second);
  for (size_t i = 0; i < view.second; i++)
    (*arr)[i] = view.first[i];

  return arr;
}
} // end extend
%enddef

// Include the real headers
%include "Teuchos_Array.hpp"
ADD_VIEW(int)
ADD_VIEW(long long)
ADD_VIEW(double)
#undef ADD_VIEW

// Make templates
%template(TeuchosArrayInt)      Teuchos::Array<int>;
%template(TeuchosArrayLongLong) Teuchos::Array<long long>;
%template(TeuchosArrayDouble)   Teuchos::Array<double>;
