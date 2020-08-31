/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include <Teuchos_ParameterList.hpp>
%}

%include <std_string.i>

// Hide warnings about overloading intrinsics
%warnfilter(314) Teuchos::ParameterList::print;
%warnfilter(314) Teuchos::ParameterList::size;

// Make the Plist an RCP
%teuchos_rcp(Teuchos::ParameterList)

// Un-camelcase the accessors
%rename(is_parameter) Teuchos::ParameterList::isParameter;
%rename(is_type) Teuchos::ParameterList::isType;

namespace Teuchos
{

class ParameterList
{
  public:
    ParameterList();
    ParameterList(const std::string& name);

    // Print the plist
    void print() const;

    template<typename T>
    void set(const std::string& name, const T& value);

    template<typename T>
    const T& get(const std::string& name);

    template<typename T>
    bool isType(const std::string& name);

    void remove(const std::string& name);

    bool isParameter(const std::string& name) const;

    ParameterList& sublist(const std::string& name);

// Instantiate get/set
%template(set) set<double>;
%template(set) set<int>;
%template(set) set<long long>;
%template(set) set<bool>;
%template(set) set<std::string>;
%template(set) set<Teuchos::Array<double> >;
%template(set) set<Teuchos::Array<int> >;
%template(set) set<Teuchos::Array<long long> >;
%template(set) set<Teuchos::ParameterList >;

%template(get_real    ) get<double>;
%template(get_integer ) get<int>;
%template(get_longlong) get<long long>;
%template(get_logical ) get<bool>;
%template(get_string  ) get<std::string>;
%template(get_arr_real    ) get<Teuchos::Array<double> >;
%template(get_arr_integer ) get<Teuchos::Array<int> >;
%template(get_arr_longlong) get<Teuchos::Array<long long> >;

}; // end class

} // end namespace Teuchos
