%module fortpetra

%import <forteuchos.i>

%{
#include "Tpetra_ConfigDefs.hpp"

typedef double                                  SC;
typedef int                                     LO;
typedef int                                     GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
%}

typedef double                                  SC;
typedef int                                     LO;
typedef int                                     GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode NO;

// Helper
namespace Kokkos {
  namespace Details {

    template<class T>
    class ArithTraits {
    public:
      typedef T val_type;
      typedef T mag_type;
    };

    template<class T>
    class InnerProductSpaceTraits {
    public:
      typedef T dot_type;
    };

  }
}
%template() Kokkos::Details::ArithTraits<SC>;
%template() Kokkos::Details::InnerProductSpaceTraits<SC>;

%ignore Teuchos::SerializationTraits;

// ignore Details namespace
%ignore Tpetra::Details;

// ignore these defines
#define TPETRA_DEPRECATED
#define KOKKOS_INLINE_FUNCTION

// Order matters!!!
%include "Tpetra_Map.i"
%include "Tpetra_Export.i"
%include "Tpetra_Import.i"
%include "Tpetra_MultiVector.i"
/* %include "Tpetra_Vector.i" */        // needs better support for inheritance
/* %include "Tpetra_Operator.i" */      // needs to understand that Tpetra::MultiVector<SC,LO,GO,NO,false> and Tpetra::MultiVector<SC,LO,GO,NO> are the same thing
%include "Tpetra_CrsGraph.i"
/* %include "Tpetra_RowMatrix.i" */     // needs better support for abstract classes
%include "Tpetra_CrsMatrix.i"
