#include "handle_helpers.hpp"
#include "fortran_operator.hpp"

namespace ForTrilinos {

  typedef double                                  SC;
  typedef int                                     LO;
  typedef int                                     GO;
  typedef Kokkos::Compat::KokkosSerialWrapperNode NO;

  auto HandleHelpers::setup_matrix_gen(const Teuchos::RCP<Teuchos::Comm<int>>& comm, int numRows, const int* rowInds, const int* rowPtrs, int numNnz, const int* colInds, const double* values) -> Teuchos::RCP<Operator> {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ArrayView;

    TEUCHOS_ASSERT(numRows >= 0);
    TEUCHOS_ASSERT((rowInds != NULL && rowPtrs != NULL) || numRows == 0);
    TEUCHOS_ASSERT(numNnz >= 0);
    TEUCHOS_ASSERT((colInds != NULL && values != NULL) || numNnz == 0);

    ArrayView<const GO> rows(rowInds, numRows);
    auto rowMap = Tpetra::createNonContigMapWithNode<LO,GO,NO>(rows, comm);

    auto A = rcp(new Matrix(rowMap, 1));

    // TODO: Can we use setAllValues?
    for (int i = 0; i < numRows; i++) {
      ArrayView<const GO> cols(colInds + rowPtrs[i], rowPtrs[i+1] - rowPtrs[i]);
      ArrayView<const SC> vals(values  + rowPtrs[i], rowPtrs[i+1] - rowPtrs[i]);

      A->insertGlobalValues(rowInds[i], cols, vals);
    }

    A->fillComplete();

    return A;
  }

  auto HandleHelpers::setup_operator_gen(const Teuchos::RCP<Teuchos::Comm<int>>& comm, int numRows, const int* rowInds, void (*funcptr)(int n, const double* x, double* y)) -> Teuchos::RCP<Operator> {
    TEUCHOS_ASSERT(numRows >= 0);
    TEUCHOS_ASSERT(rowInds != NULL || numRows == 0);

    // NOTE: we make a major assumption on the maps:
    //      rowMap == domainMap == rangeMap
    Teuchos::ArrayView<const GO> rows(rowInds, numRows);
    auto map = Tpetra::createNonContigMapWithNode<LO,GO,NO>(rows, comm);

    return Teuchos::rcp(new FortranOperator(funcptr, map, map));
  }

}
