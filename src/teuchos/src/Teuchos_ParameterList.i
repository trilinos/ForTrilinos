/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include "Teuchos_ParameterList.hpp"
%}

%include <std_string.i>

// Instantiate std::pair view
%include <typemaps.i>
%fortran_string_view(char)
%fortran_view(int)
%fortran_view(double)

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

// Construct without a name
ParameterList()
{
    return new Teuchos::ParameterList();
}

// Constructor
ParameterList(std::pair<const char*, size_t> name)
{
    return new Teuchos::ParameterList(std::string(name.first, name.second));
}

// >>> SCALAR SET/GET

template<typename T>
void set_scalar(std::pair<const char*, size_t> name, const T& value)
{
    $self->set(std::string(name.first, name.second), value);
}

template<typename T>
void get_scalar(std::pair<const char*, size_t> name, T& value)
{
    value = $self->get<T>(std::string(name.first, name.second));
}

// Instantiate get/set
%template(get) get_scalar<double>;
%template(set) set_scalar<double>;
%template(get) get_scalar<int>;
%template(set) set_scalar<int>;

// >>> STRING SET/GET

void set(std::pair<const char*, size_t> name,
         std::pair<const char*, size_t> value)
{
    $self->set(std::string(name.first,  name.second),
               std::string(value.first, value.second));
}

void get(std::pair<const char*, size_t> name,
         std::pair<char*, size_t>       value)
{
    const std::string& str
        = $self->get<std::string>(std::string(name.first, name.second));

    // FIXME: should we check the size here?
    if (value.second < (size_t)str.size()) {
        std::ostringstream os;
        os << "Array size mismatch: " << str.size() << " > " << value.second;
        throw std::range_error(os.str());
    }
    std::copy(str.begin(), str.end(), value.first);
    std::fill_n(value.first+str.size(), value.second - str.size(), ' ');
}

// >>> ARRAY SET/GET
template<typename T>
void set_array(std::pair<const char*, size_t> name,
               std::pair<const T*, size_t>    value)
{
    typedef Teuchos::Array<T> ArrayT;
    $self->set(std::string(name.first, name.second),
               ArrayT(value.first, value.first + value.second));
}

template<typename T>
void get_array(std::pair<const char*, size_t> name,
               std::pair<const T*, size_t>    value)
{
    typedef Teuchos::Array<T> ArrayT;
    const ArrayT& arr
        = $self->get<ArrayT>(std::string(name.first, name.second));

    // FIXME: should we check the size here?
    if (value.second < (size_t)arr.size()) {
        std::ostringstream os;
        os << "Array size mismatch: " << arr.size() << " != " << value.second;
        throw std::range_error(os.str());
    }


    value.first = arr.getRawPtr();
    value.second = arr.size();
}

%template(set) set_array<double>;
%template(get) get_array<double>;
%template(set) set_array<int>;
%template(get) get_array<int>;

// >>> PLIST SET/GET

void set(std::pair<const char*, size_t> name,
         Teuchos::RCP<Teuchos::ParameterList> plist)
{
    $self->set(std::string(name.first, name.second), *plist);
}

void get(std::pair<const char*, size_t> name,
         Teuchos::RCP<Teuchos::ParameterList>& plist)
{
    plist = Teuchos::sublist(
            Teuchos::rcpFromRef(*$self),
            std::string(name.first, name.second),
            true); // must exist
}

Teuchos::RCP<Teuchos::ParameterList> sublist(std::pair<const char*, size_t> name)
{
    return Teuchos::sublist(
            Teuchos::rcpFromRef(*$self),
            std::string(name.first, name.second));
}

// >>> TYPE-FREE QUERIES

//! Get the length of a string or array
int get_length(std::pair<const char*, size_t> name)
{
    std::string key(name.first, name.second);

#define PLIST_CHECK_RETURN(TYPE) \
  if ($self->isType<TYPE >(key)) \
    return $self->get<TYPE >(key).size();

    PLIST_CHECK_RETURN(std::string);
    PLIST_CHECK_RETURN(Teuchos::Array<int>)
    PLIST_CHECK_RETURN(Teuchos::Array<double>)
#undef PLIST_CHECK_RETURN

    // No type found
    return -1;
}

//! Delete a parameter
void remove(std::pair<const char*, size_t> name)
{
    $self->remove(std::string(name.first, name.second));
}

//! Whether the given parameter exists
bool is_parameter(std::pair<const char*, size_t> name) const
{
    return $self->isParameter(std::string(name.first, name.second));
}

} // End %extend
}; // end class

} // end namespace Teuchos
