/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include "Teuchos_ArrayRCP.hpp"
%}

%include <view.i>

namespace Teuchos
{
template<typename T>
class ArrayRCP
{
  public:
    typedef T value_type;
    typedef std::size_t size_type;
};
}

// Create RCP with no-ownership flag (it's set to false)
%define TEUCHOS_ARRAYRCP_IMPL(CONST, CTYPE)
  // Treat const references as values
  %apply Teuchos::ArrayRCP<CONST CTYPE> { const Teuchos::ArrayRCP<CONST CTYPE>& }

  // C input translation typemaps: $1 is ArrayRCP, $input is SwigArrayWrapper
  %typemap(in) Teuchos::ArrayRCP<CONST CTYPE> %{
    $1 = Teuchos::ArrayRCP(static_cast<CONST CTYPE*>($input->data), 0, $input->size, false);
  %}
  %typemap(in) const Teuchos::ArrayRCP<CONST CTYPE>& (Teuchos::ArrayRCP<CONST CTYPE> temprcp) %{
    temprcp = Teuchos::ArrayRCP(static_cast<CONST CTYPE*>($input->data), 0, $input->size, false);
    $1 = &temprcp;
  %}

  %template() Teuchos::ArrayRCP<CONST CTYPE>;
%enddef

/* -------------------------------------------------------------------------
 * \def %TEUCHOS_ARRAYRCP
 *
 * This maps a return value of ArrayRCP to a small struct (mirrored in
 * fortran) that defines the start and size of a contiguous array.
 */
%define TEUCHOS_ARRAYRCP(CTYPE)
  // Use SwigArrayWrapper and array pointer translation for this type
  FORT_ARRAYPTR_TYPEMAP(CTYPE, Teuchos::ArrayRCP<CTYPE>)
  %apply Teuchos::ArrayRCP<CTYPE> { Teuchos::ArrayRCP<const CTYPE> }

  // Use pointer translation for const variables too
  TEUCHOS_ARRAYRCP_IMPL(, CTYPE)
  TEUCHOS_ARRAYRCP_IMPL(const, CTYPE)

  %typemap(out) Teuchos::ArrayRCP<const CTYPE> %{
    $result.data = const_cast<CTYPE*>($1.getRawPtr());
    $result.size = $1.size();
  %}
%enddef

TEUCHOS_ARRAYRCP(int);
TEUCHOS_ARRAYRCP(double);
TEUCHOS_ARRAYRCP(long long);
TEUCHOS_ARRAYRCP(size_t);

