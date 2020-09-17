/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include <Teuchos_ArrayView.hpp>
%}

%include <fortranarray.swg>

namespace Teuchos
{
template<typename _Tp>
class ArrayView
{
    public:
    typedef _Tp value_type;
    typedef _Tp* pointer;
    typedef std::size_t size_type;

    // Add native wrapping typemaps to convert to/from Teuchos array
    %fortran_array_pointer(_Tp, ArrayView<_Tp>)

    %typemap(in, noblock=1) ArrayView<_Tp> {
      $1 = Teuchos::ArrayView<_Tp>(static_cast<_Tp*>($input->data), $input->size);
    }

    %typemap(out, noblock=1) ArrayView<_Tp> {
      $result.data = (void*)$1.getRawPtr();
      $result.size = $1.size();
    }

    %apply ArrayView<_Tp> { const ArrayView<_Tp> & }

    %typemap(in, noblock=1) const ArrayView<_Tp> & (Teuchos::ArrayView<_Tp> tmpview) {
      tmpview = Teuchos::ArrayView<_Tp>(static_cast<_Tp*>($input->data), $input->size);
      $1 = &tmpview;
    }
  %typemap(out, noblock=1) const ArrayView<_Tp>& {
    $result.data = (void*)($1->getRawPtr());
    $result.size = $1->size();
  }

    // ...and by reference, modifying an incoming pointer
    %fortran_array_handle(_Tp, ArrayView<_Tp>)
  %typemap(in, noblock=1) ArrayView<_Tp>& (ArrayView<_Tp> tmpview) {
    tmpview = $1_basetype(static_cast<$1_basetype::pointer>($input->data), $input->size);
    $1 = &tmpview;
  }
  %typemap(argout, noblock=1, match="in") ArrayView<_Tp>& {
    $input->data = (void*)tmpview$argnum.getRawPtr();
    $input->size = tmpview$argnum.size();
  }
};
}

// Instantiate typemaps for both TYPE and const TYPE
%define %teuchos_instantiate_arrayview_typemaps(TYPE)
%template() Teuchos::ArrayView<TYPE>;
%template() Teuchos::ArrayView<const TYPE>;
%enddef


%teuchos_instantiate_arrayview_typemaps(int);
%teuchos_instantiate_arrayview_typemaps(double);
%teuchos_instantiate_arrayview_typemaps(long long);
%teuchos_instantiate_arrayview_typemaps(unsigned int);
%teuchos_instantiate_arrayview_typemaps(unsigned long);
%teuchos_instantiate_arrayview_typemaps(unsigned long long);
