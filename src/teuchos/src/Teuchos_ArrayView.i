// Dependencies
%include "Teuchos_RCP.i"

%{
#include "Teuchos_ArrayViewDecl.hpp"
%}

// Typedefs
typedef int Teuchos_Ordinal;

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

// Include the real headers
%include "Teuchos_ArrayViewDecl.hpp"

// Make templates
%template(TeuchosArrayViewInt)          Teuchos::ArrayView<int>;
%template(TeuchosArrayViewIntConst)     Teuchos::ArrayView<const int>;
%template(TeuchosArrayViewDouble)       Teuchos::ArrayView<double>;
%template(TeuchosArrayViewDoubleConst)  Teuchos::ArrayView<const double>;
%template(TeuchosArrayViewSizeT)        Teuchos::ArrayView<size_t>;
%template(TeuchosArrayViewSizeTConst)   Teuchos::ArrayView<const size_t>;
