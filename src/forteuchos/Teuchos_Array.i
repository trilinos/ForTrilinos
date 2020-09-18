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
%include <std_common.i>
%include <typemaps.i>

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

    %extend {
  %apply (const SWIGTYPE *DATA, ::size_t SIZE)
    { (const _Tp* DATA, size_type SIZE) };

  // Assign from another Array
  void assign(const _Tp* DATA, size_type SIZE) {
    $self->assign(DATA, DATA + SIZE);
  }

  // Convert views to and from native Fortran arrays/pointers
  %fortran_array_pointer(_Tp, Array<_Tp>& view);

  %typemap(in, noblock=1) Array<_Tp>& view (Array<_Tp> temp){
    temp = $*1_ltype(static_cast<$1_basetype::pointer>($input->data),
                   static_cast<$1_basetype::pointer>($input->data)
                   + $input->size);
    $1 = &temp;
  }

  %typemap(out, noblock=1) Array<_Tp>& view {
    $result.data = ($1->empty() ? NULL : &(*$1->begin()));
    $result.size = $1->size();
  }

  Array<_Tp>& view() {
    return *$self;
  }
    }
};
}

%template(TeuchosArrayInt) Teuchos::Array<int>;
%template(TeuchosArrayDbl) Teuchos::Array<double>;
%template(TeuchosArrayLongLong) Teuchos::Array<long long>;
