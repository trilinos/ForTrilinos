/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include "Teuchos_ArrayRCP.hpp"
%}

%include <fortranarray.swg>

namespace Teuchos
{
template<typename _Tp>
class ArrayRCP
{
  public:
    typedef _Tp value_type;
    typedef std::size_t size_type;

#if 0
    // Rename 'call' operator to 'view'
    %rename(view) operator();

    // Add accessors to get a view
    ArrayRCP<_Tp> operator()(size_type lowerOffset, size_type size) const;
    ArrayRCP<_Tp> operator()() const;

    // Smart pointer attributes
    int strong_count() const;
    int weak_count() const;
    int total_count() const;
    bool has_ownership() const;
#endif

    // Add native wrapping typemaps to convert *input* to Teuchos RCP
    %fortran_array_pointer(_Tp, ArrayRCP<_Tp>)
    %apply ArrayRCP<_Tp> { const ArrayRCP<_Tp> & }

    // Output isn't currently implemented
    %typemap(out) ArrayRCP<_Tp>
    "#error Can't return $type by value"
    %typemap(out) const ArrayRCP<_Tp>&
    "#error Can't return $type by const reference"

    // Wrap RCP by value
    %typemap(in, noblock=1) ArrayRCP<_Tp> {
      $1 = Teuchos::ArrayRCP<_Tp>(static_cast<_Tp*>($input->data), 0, $input->size, false, Teuchos::RCP_DISABLE_NODE_LOOKUP);
    }
    // Wrap RCP by const reference
    %typemap(in, noblock=1) const ArrayRCP<_Tp>& (ArrayRCP<_Tp> tmparr) {
      tmparr = Teuchos::ArrayRCP<_Tp>(static_cast<_Tp*>($input->data), 0, $input->size, false, Teuchos::RCP_DISABLE_NODE_LOOKUP);
      $1 = &tmparr;
    }

};
}

%template() Teuchos::ArrayRCP<int>;
%template() Teuchos::ArrayRCP<double>;
%template() Teuchos::ArrayRCP<long long>;

%template() Teuchos::ArrayRCP<unsigned int>;
%template() Teuchos::ArrayRCP<unsigned long>;
%template() Teuchos::ArrayRCP<unsigned long long>;

%template() Teuchos::ArrayRCP<const int>;
%template() Teuchos::ArrayRCP<const double>;
%template() Teuchos::ArrayRCP<const long long>;

%template() Teuchos::ArrayRCP<const unsigned int>;
%template() Teuchos::ArrayRCP<const unsigned long>;
%template() Teuchos::ArrayRCP<const unsigned long long>;

#if 0
%define %teuchos_template_arrayrcp(EXT, CPPTYPE)

// Only instantiate non-const
%template(ArrayRCP ## EXT) Teuchos::ArrayRCP<CPPTYPE>;

// Use non-const typemaps for RCP-to-const
%apply Teuchos::ArrayRCP<CPPTYPE> { Teuchos::ArrayRCP<const CPPTYPE> };

// Use arcp const cast to convert mutable to const
%typemap(in, noblock=1) Teuchos::ArrayRCP<const CPPTYPE>
(Teuchos::ArrayRCP<CPPTYPE>* origrcp) {
  SWIG_check_nonnull(*$input, "$1_ltype", "$fclassname", "$decl", return $null);
  origrcp = %static_cast($input->cptr, Teuchos::ArrayRCP<const CPPTYPE>*);
  $1 = Teuchos::arcp_const_cast<const CPPTYPE>(*origrcp);
}

%enddef

// Instantiate typemaps only
%template(ArrayRCP_Int) Teuchos::ArrayRCP<int>;
%template(ArrayRCP_Dbl) Teuchos::ArrayRCP<double>;
%template(ArrayRCP_Longlong) Teuchos::ArrayRCP<long long>;
#endif

