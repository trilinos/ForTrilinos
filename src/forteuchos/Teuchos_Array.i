/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include <Teuchos_Array.hpp>
%}

// Include typemap for const array references
%include <std_container.i>

namespace Teuchos
{
template<typename _Tp>
class Array
{
  public:
    typedef _Tp value_type;
    typedef std::size_t size_type;

    // Implicit conversion from array view
    Array(const ArrayView<const _Tp>&);

    // Wrap like a std::vector: const references and values are treated as
    // native Fortran data; mutable references and pointers are wrapped.
    %std_native_container(Teuchos::Array<_Tp>)

    %extend {
        ArrayView<_Tp> view() {
            return (*$self)();
        }
    }
};
}

%template(TeuchosArrayInt) Teuchos::Array<int>;
%template(TeuchosArrayDbl) Teuchos::Array<double>;
%template(TeuchosArrayLongLong) Teuchos::Array<long long>;
