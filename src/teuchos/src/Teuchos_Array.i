%module forteuchos_array

// Dependencies
%include "Teuchos_RCP.i"

%{
#include "Teuchos_Array.hpp"
%}

// Typedefs
typedef int Teuchos_Ordinal;

// Make RCPs
%teuchos_rcp(Teuchos::Array<int>)
%teuchos_rcp(Teuchos::Array<double>)

// Ignore misc classes
%ignore TypeNameTraits;
%ignore InvalidArrayStringRepresentation;

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

// Include the real headers
%include "Teuchos_Array.hpp"

// Make templates
%template(TeuchosArrayInt)      Teuchos::Array<int>;
%template(TeuchosArrayDouble)   Teuchos::Array<double>;
