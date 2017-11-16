// Dependencies
%include "Teuchos_RCP.i"

%{
#include "Teuchos_ArrayViewDecl.hpp"
%}

%include <typemaps.i>

// Make RCPs
%teuchos_rcp(Teuchos::ArrayView<int>)
%teuchos_rcp(Teuchos::ArrayView<const int>)
%teuchos_rcp(Teuchos::ArrayView<double>)
%teuchos_rcp(Teuchos::ArrayView<const double>)
%teuchos_rcp(Teuchos::ArrayView<size_t>)
%teuchos_rcp(Teuchos::ArrayView<const size_t>)

// Ignore misc classes
/* %ignore TypeNameTraits; */

// Ignore functions
%ignore Teuchos::ArrayView::assign;
%ignore Teuchos::ArrayView::begin;
%ignore Teuchos::ArrayView::end;
%ignore Teuchos::ArrayView::ArrayView( ENull null_arg = null );
%ignore Teuchos::ArrayView::ArrayView(T* p, size_type size, const ERCPNodeLookup rcpNodeLookup); // ignore the non-default ctor
%ignore Teuchos::ArrayView::ArrayView(const T* p, size_type size, const ERCPNodeLookup rcpNodeLookup); // ignore the non-default ctor
%ignore Teuchos::ArrayView::access_private_ptr;
%ignore Teuchos::ArrayView::access_private_arcp;
%ignore Teuchos::ArrayView::getRawPtr;

// Postpone
%ignore Teuchos::ArrayView(std::vector<typename ConstTypeTraits<T>::NonConstType>& vec);
%ignore Teuchos::ArrayView::operator[];
%ignore Teuchos::ArrayView::operator();
%ignore Teuchos::ArrayView::operator=(const ArrayView<T>& array);
%ignore Teuchos::ArrayView::operator=(const ArrayView<const T>& array);
%ignore Teuchos::ArrayView::ArrayView(std::vector<typename ConstTypeTraits<T>::NonConstType>& vec);
%ignore Teuchos::ArrayView::ArrayView(const std::vector<typename ConstTypeTraits<T>::NonConstType>& vec);
%ignore Teuchos::ArrayView::getConst() const;
%ignore Teuchos::ArrayView::operator ArrayView<const T>() const;

%define ADD_VIEW(TYPE)
// Extend ArrayView
%extend Teuchos::ArrayView<TYPE> {
ArrayView(std::pair<TYPE*, size_t> view)
{
  return new Teuchos::ArrayView<TYPE>(view.first, view.second);
}

  /* std::pair<TYPE*, std::size_t> view() { */
    /* if ($self->size()) */
      /* return {nullptr, 0}; */
    /* return {$self->getRawPtr(), $self->size()}; */
  /* } */
} // end extend
%enddef

// Include the real headers
%include "Teuchos_ArrayViewDecl.hpp"
ADD_VIEW(int)
ADD_VIEW(const int)
ADD_VIEW(long long)
ADD_VIEW(const long long)
ADD_VIEW(double)
ADD_VIEW(const double)
ADD_VIEW(size_t)
ADD_VIEW(const size_t)
#undef ADD_VIEW

// Make templates
%template(TeuchosArrayViewInt)          Teuchos::ArrayView<int>;
%template(TeuchosArrayViewIntConst)     Teuchos::ArrayView<const int>;
%template(TeuchosArrayViewLongLong)     Teuchos::ArrayView<long long>;
%template(TeuchosArrayViewLongLongConst)Teuchos::ArrayView<const long long>;
%template(TeuchosArrayViewDouble)       Teuchos::ArrayView<double>;
%template(TeuchosArrayViewDoubleConst)  Teuchos::ArrayView<const double>;
%template(TeuchosArrayViewSizeT)        Teuchos::ArrayView<size_t>;
%template(TeuchosArrayViewSizeTConst)   Teuchos::ArrayView<const size_t>;
