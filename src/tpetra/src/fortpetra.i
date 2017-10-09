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
  }
}
%template() Kokkos::Details::ArithTraits<SC>;

// ignore Details namespace
%ignore Tpetra::Details;

// Order matters!!!
%include "Tpetra_Map.i"
%include "Tpetra_MultiVector.i"
%include "Tpetra_CrsMatrix.i"