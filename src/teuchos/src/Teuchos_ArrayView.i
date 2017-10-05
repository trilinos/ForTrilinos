// Dependencies
%include <std_vector.i>
%include "Teuchos_RCP.i"

%{
#include "Teuchos_ArrayViewDecl.hpp"
%}

// Typedefs
typedef int Teuchos_Ordinal;

%template(VectorInt)    std::vector<int>;
%template(VectorDouble) std::vector<double>;

// Make RCPs
%teuchos_rcp(Teuchos::ArrayView<int>)
%teuchos_rcp(Teuchos::ArrayView<double>)

// Ignore misc classes
/* %ignore TypeNameTraits; */

// Ignore functions
/* %ignore Teuchos::ArrayView::assign; */
%ignore Teuchos::ArrayView::assign;
%ignore Teuchos::ArrayView::begin;
%ignore Teuchos::ArrayView::end;
%ignore Teuchos::ArrayView::ArrayView( ENull null_arg = null );
%ignore Teuchos::ArrayView::ArrayView(T* p, size_type size, const ERCPNodeLookup rcpNodeLookup); // ignore the non-default ctor

// Postpone
/* %ignore Teuchos::ArrayView::operator=(const Array<T>& a); */
%ignore Teuchos::ArrayView(std::vector<typename ConstTypeTraits<T>::NonConstType>& vec);
%ignore Teuchos::ArrayView::operator[];
%ignore Teuchos::ArrayView::operator();
%ignore Teuchos::ArrayView::operator=(const ArrayView<T>& array);
%ignore Teuchos::ArrayView::ArrayView(std::vector<typename ConstTypeTraits<T>::NonConstType>& vec);
%ignore Teuchos::ArrayView::ArrayView(const std::vector<typename ConstTypeTraits<T>::NonConstType>& vec);
%ignore Teuchos::ArrayView::getConst() const;
%ignore Teuchos::ArrayView::operator ArrayView<const T>() const;

// Include the real headers
%include "Teuchos_ArrayViewDecl.hpp"

// Make templates
%template(TeuchosArrayViewInt)      Teuchos::ArrayView<int>;
