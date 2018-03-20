/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include "Teuchos_ArrayView.hpp"
%}

%include <view.i>

namespace Teuchos
{
template<typename T>
class ArrayView
{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
};
}

%define TEUCHOS_ARRAYVIEW_IMPL(CONST, CTYPE)
  // Treat const references as values
  %apply Teuchos::ArrayView<CONST CTYPE> { const Teuchos::ArrayView<CONST CTYPE>& }

  // C input translation typemaps: $1 is ArrayView, $input is SwigArrayWrapper
  %typemap(in) Teuchos::ArrayView<CONST CTYPE> %{
    $1 = Teuchos::arrayView(static_cast<CONST CTYPE*>($input->data), $input->size);
  %}
  %typemap(in) const Teuchos::ArrayView<CONST CTYPE>& (Teuchos::ArrayView<CONST CTYPE> tempview) %{
    tempview = Teuchos::arrayView(static_cast<CONST CTYPE*>($input->data), $input->size);
    $1 = &tempview;
  %}

  %template() Teuchos::ArrayView<CONST CTYPE>;
%enddef

/* -------------------------------------------------------------------------
 * \def %TEUCHOS_ARRAYVIEW
 *
 * This maps a return value of ArrayView to a small struct (mirrored in
 * fortran) that defines the start and size of a contiguous array.
 */
%define TEUCHOS_ARRAYVIEW(CTYPE)
  // Use SwigArrayWrapper and array pointer translation for this type
  FORT_ARRAYPTR_TYPEMAP(CTYPE, Teuchos::ArrayView<CTYPE>)
  %apply Teuchos::ArrayView<CTYPE> { Teuchos::ArrayView<const CTYPE> }

  // Use pointer translation for const variables too
  TEUCHOS_ARRAYVIEW_IMPL(, CTYPE)
  TEUCHOS_ARRAYVIEW_IMPL(const, CTYPE)

  %typemap(out) Teuchos::ArrayView<const CTYPE> %{
    $result.data = const_cast<CTYPE*>($1.getRawPtr());
    $result.size = $1.size();
  %}
%enddef

TEUCHOS_ARRAYVIEW(int);
TEUCHOS_ARRAYVIEW(double);
TEUCHOS_ARRAYVIEW(long long);
TEUCHOS_ARRAYVIEW(size_t);


/* -------------------------------------------------------------------------
 * Typemaps: convert to/from Fortran indexing.
 */
%typemap(in, noblock=1) Teuchos::ArrayView<int> FORTRAN_INDEX
    (Teuchos::Array<$1_ltype::value_type > tmparr,
     Teuchos::ArrayView<$1_ltype::value_type > tmpview) {
  tmparr.resize($1->size());

  for (int i = 0; i < $1->size(); i++)
    (*$1)[i] = (static_cast<$1_ltype::value_type*>($input))[i] - 1;
  tmpview = tmparr
  $1 = &tmpview;
}

%typemap(argout, noblock=1) Teuchos::ArrayView<int> FORTRAN_INDEX {
  for (int i = 0; i < $1->size(); i++)
    (*$1)[i]++;
}

// Convert a return-by-view array in C indexing to Fortran indexing.
%typemap(argout, noblock=1) Teuchos::ArrayView<int> FORTRAN_INDEX {
  for (int i = 0; i < $1->size(); i++)
    (*$1)[i]++;
}
