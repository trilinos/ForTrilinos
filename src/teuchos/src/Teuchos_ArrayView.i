/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include "Teuchos_ArrayView.hpp"
%}

%include <fortranarray.swg>

namespace Teuchos
{
template<typename _Tp>
class ArrayView
{
    public:
    typedef _Tp value_type;
    typedef std::size_t size_type;

    // Add native wrapping typemaps to convert to/from Teuchos array
    %fortran_array_pointer(_Tp, ArrayView<_Tp>)

    %typemap(in, noblock=1) ArrayView<_Tp> (_Tp* tempbegin) {
      tempbegin = static_cast<_Tp*>($input->data);
      $1 = Teuchos::ArrayView<_Tp>(tempbegin, $input->size);
    }

    %typemap(in, noblock=1) const ArrayView<_Tp> & (_Tp* tempbegin, Teuchos::ArrayView<_Tp> temparr) {
      tempbegin = static_cast<_Tp*>($input->data);
      temparr = Teuchos::ArrayView<_Tp>(tempbegin, $input->size);
      $1 = &temparr;
    }

    %typemap(out, noblock=1) ArrayView<_Tp> {
      $result.data = (void*)$1.getRawPtr();
      $result.size = $1.size();
    }

    // XXX: %apply only works for *NAME-INSTANTIATED* types, not %template()
    // types, so we have to manually copy these typemaps.
    %typemap(ctype) const ArrayView<_Tp> & = ArrayView<_Tp>;
    %typemap(imtype) const ArrayView<_Tp> & = ArrayView<_Tp>;
    %typemap(ftype) const ArrayView<_Tp> & = ArrayView<_Tp>;
    %typemap(fin) const ArrayView<_Tp> & = ArrayView<_Tp>;
    %typemap(findecl) const ArrayView<_Tp> & = ArrayView<_Tp>;
};
}

// Instantiate typemaps only
%template() Teuchos::ArrayView<int>;
%template() Teuchos::ArrayView<double>;
%template() Teuchos::ArrayView<long long>;

%template() Teuchos::ArrayView<const int>;
%template() Teuchos::ArrayView<const double>;
%template() Teuchos::ArrayView<const long long>;

%template() Teuchos::ArrayView<const unsigned int>;
%template() Teuchos::ArrayView<const unsigned long>;
%template() Teuchos::ArrayView<const unsigned long long>;
