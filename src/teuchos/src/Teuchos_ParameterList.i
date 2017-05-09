//---------------------------------*-SWIG-*----------------------------------//
/*!
 * \file   parameterlist/Teuchos_ParameterList.i
 * \author Seth R Johnson
 * \date   Tue Dec 06 17:59:15 2016
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
%{
#include "Teuchos_ParameterList.hpp"
%}

%include <std_string.i>
%include <typemaps.i>

// Use string copying and size checking
%fragment("StringCopyout");
%fragment("ArraySizeCheck");

// Hide warnings about overloading intrinsics
%warnfilter(314) Teuchos::ParameterList::print;
%warnfilter(314) Teuchos::ParameterList::size;

// Make the Plist an RCP
%teuchos_rcp(Teuchos::ParameterList)


namespace Teuchos
{

class ParameterList
{
  public:

    // Print the plist
    void print() const;

    %extend {

// Constructor
ParameterList(const char* STRING, int SIZE)
{
    return new Teuchos::ParameterList(std::string(STRING, SIZE));
}

// >>> TYPED QUERIES

// Use fortran string typemap for string value as well as key

//! String set/get
template<typename T>
void set_scalar(const char* STRING, int SIZE, const T& value)
{
    $self->set(std::string(STRING, SIZE), value);
}

template<typename T>
void get_scalar(const char* STRING, int SIZE, T& value)
{
    value = $self->get<T>(std::string(STRING, SIZE));
}

// Instantiate get/set
%template(get) get_scalar<double>;
%template(set) set_scalar<double>;
%template(get) get_scalar<int>;
%template(set) set_scalar<int>;
%template(get) get_scalar<Teuchos::ParameterList>;
%template(set) set_scalar<Teuchos::ParameterList>;

%apply (char* STRING, int SIZE)
{ (char* VALSTRING, int VALSIZE) };
%apply (const char* STRING, int SIZE)
{ (const char* VALSTRING, int VALSIZE) };

// String get/set
void set(const char* STRING, int SIZE,
         const char* VALSTRING, int VALSIZE)
{
    $self->set(std::string(STRING, SIZE),
               std::string(VALSTRING, VALSIZE));
}

void get(const char* STRING, int SIZE,
         char* VALSTRING,    int VALSIZE)
{
    const std::string& value
        = $self->get<std::string>(std::string(STRING, SIZE));
    swig::string_copyout(value, VALSTRING, VALSIZE);
}

// Use fortran array typemap for value as well as key
%apply (SWIGTYPE* ARRAY, int SIZE)
{ (int* ARRAY, int ARRAYSIZE),
  (double* ARRAY, int ARRAYSIZE),
  (const int* ARRAY, int ARRAYSIZE),
  (const double* ARRAY, int ARRAYSIZE) };

//! Array set
template<typename T>
void set_array(const char* STRING, int SIZE,
               const T* ARRAY,     int ARRAYSIZE)
{
    typedef Teuchos::Array<T> ArrayT;
    $self->set(std::string(STRING, SIZE),
               ArrayT(ARRAY, ARRAY + ARRAYSIZE));
}

//! Array get
template<typename T>
void get_array(const char* STRING, int SIZE,
               T* ARRAY,           int ARRAYSIZE)
{
    typedef Teuchos::Array<T> ArrayT;
    const ArrayT& arr
        = $self->get<ArrayT>(std::string(STRING, SIZE));

    swig::array_size_check(arr.size(), ARRAYSIZE);
    std::copy(arr.begin(), arr.end(), ARRAY);
}

%template(set) set_array<double>;
%template(get) get_array<double>;
%template(set) set_array<int>;
%template(get) get_array<int>;

//! Get the length of a string or array
int get_length(const char* STRING, int SIZE)
{
    std::string key(STRING, SIZE);
#define PLIST_CHECK_RETURN(TYPE) \
if ($self->isType<TYPE >(key)) return $self->get<TYPE >(key).size();

    PLIST_CHECK_RETURN(std::string);
    PLIST_CHECK_RETURN(Teuchos::Array<int>)
    PLIST_CHECK_RETURN(Teuchos::Array<double>)

#undef PLIST_CHECK_RETURN
    // No type found
    return -1;
}

//! Delete a parameter
void remove(const char* STRING, int SIZE)
{
    $self->remove(std::string(STRING, SIZE));
}

//! Whether the given parameter exists
bool is_parameter(const char* STRING, int SIZE) const
{
    return $self->isParameter(std::string(STRING, SIZE));
}


} // End %extend
}; // end class

} // end namespace Teuchos

//---------------------------------------------------------------------------//
// end of parameterlist/Teuchos_ParameterList.i
//---------------------------------------------------------------------------//
