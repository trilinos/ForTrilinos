#ifndef FORTRILINOS_HANDLE_HELPERS_HPP
#define FORTRILINOS_HANDLE_HELPERS_HPP

#include <Tpetra_Operator.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_RCP.hpp>

#include "fortran_operator.hpp"

namespace ForTrilinos {


  class HandleHelpers {
  public:
    typedef double                                  SC;
    typedef int                                     LO;
    typedef long long                               GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;

  public:
    typedef Tpetra::Operator<SC,LO,GO,NO>           Operator;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO>          Matrix;

  public:

    // Setup matrix
    static Teuchos::RCP<Matrix>
    setup_matrix_gen(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                     std::pair<const GO*,size_t> rowInds, std::pair<const LO*,size_t> rowPtrs,
                     std::pair<const GO*,size_t> colInds, std::pair<const SC*,size_t> values);

    // Setup operator
    static Teuchos::RCP<Operator>
    setup_operator_gen(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                       std::pair<const GO*, size_t> rowInds, OperatorCallback callback);
  };

}

#endif // FORTRILINOS_HANDLE_HELPERS_HPP
