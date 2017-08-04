#ifndef FORTRILINOS_HANDLE_HELPERS_HPP
#define FORTRILINOS_HANDLE_HELPERS_HPP

#include <Tpetra_Operator.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_RCP.hpp>

namespace ForTrilinos {


  class HandleHelpers {
  private:
    typedef double                                  SC;
    typedef int                                     LO;
    typedef int                                     GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;

  public:
    typedef Tpetra::Operator<SC,LO,GO,NO>           Operator;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO>          Matrix;
    typedef void (*OperatorCallback)(int n, const double* x, double* y);

  public:

    // Setup matrix
    static Teuchos::RCP<Operator>
    setup_matrix_gen(const Teuchos::RCP<Teuchos::Comm<int>>& comm, int numRows, const int* rowInds, const int* rowPtrs, int numNnz, const int* colInds, const double* values);

    // Setup operator
    static Teuchos::RCP<Operator>
    setup_operator_gen(const Teuchos::RCP<Teuchos::Comm<int>>& comm, int numRows, const int* rowInds, OperatorCallback callback);
  };

}

#endif // FORTRILINOS_HANDLE_HELPERS_HPP
