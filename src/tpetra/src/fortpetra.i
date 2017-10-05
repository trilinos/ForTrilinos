%module fortpetra

%import <forteuchos.i>

%{
#include "Tpetra_ConfigDefs.hpp"

typedef double Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;
typedef Kokkos::Compat::KokkosSerialWrapperNode Node;
%}

typedef double Scalar;
typedef int LocalOrdinal;
typedef int GlobalOrdinal;
typedef Kokkos::Compat::KokkosSerialWrapperNode Node;

// Helper
namespace Kokkos {
  namespace Details {
    template<class T>
    class ArithTraits {
    public:
      typedef T val_type;
      typedef T mag_type;
    };
  }
}
%template(ArithTraitsScalar) Kokkos::Details::ArithTraits<Scalar>;

// Order matters!!!
%include "Tpetra_Map.i"
%include "Tpetra_MultiVector.i"
%include "Tpetra_CrsMatrix.i"
