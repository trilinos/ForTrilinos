#include "handle_helpers.hpp"
#include "fortran_operator.hpp"

namespace ForTrilinos {

  auto HandleHelpers::setup_matrix_gen(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
        std::pair<const GO*,size_t> rowInds_pair, std::pair<const LO*,size_t> rowPtrs_pair,
        std::pair<const GO*,size_t> colInds_pair, std::pair<const SC*,size_t> values_pair) -> Teuchos::RCP<Matrix> {
    using Teuchos::ArrayView;

    auto numRows = rowInds_pair.second;
    auto numNnz  = colInds_pair.second;

    TEUCHOS_ASSERT(rowPtrs_pair.second == numRows+1);
    TEUCHOS_ASSERT(values_pair .second == numNnz);

    auto rowInds = rowInds_pair.first;
    auto rowPtrs = rowPtrs_pair.first;
    auto colInds = colInds_pair.first;
    auto values  = values_pair .first;

    TEUCHOS_ASSERT(numRows >= 0);
    TEUCHOS_ASSERT(numNnz  >= 0);
    TEUCHOS_ASSERT((rowInds != NULL && rowPtrs != NULL) || numRows == 0);
    TEUCHOS_ASSERT((colInds != NULL && values != NULL)  || numNnz  == 0);

    ArrayView<const GO> rows(rowInds, numRows);
    auto rowMap = Tpetra::createNonContigMapWithNode<LO,GO,NO>(rows, comm);

    auto A = Teuchos::rcp(new Matrix(rowMap, 1));

    // TODO: Can we use setAllValues?
    for (int i = 0; i < numRows; i++) {
      ArrayView<const GO> cols(colInds + rowPtrs[i], rowPtrs[i+1] - rowPtrs[i]);
      ArrayView<const SC> vals(values  + rowPtrs[i], rowPtrs[i+1] - rowPtrs[i]);

      A->insertGlobalValues(rowInds[i], cols, vals);
    }

    A->fillComplete();

    return A;
  }

  auto HandleHelpers::setup_operator_gen(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
        std::pair<const GO*, size_t> rowInds_pair, OperatorCallback callback) -> Teuchos::RCP<Operator> {
    auto rowInds = rowInds_pair.first;
    auto numRows = rowInds_pair.second;
    TEUCHOS_ASSERT(numRows >= 0);
    TEUCHOS_ASSERT(rowInds != NULL || numRows == 0);

    // NOTE: we make a major assumption on the maps:
    //      rowMap == domainMap == rangeMap
    Teuchos::ArrayView<const GO> rows(rowInds, numRows);
    auto map = Tpetra::createNonContigMapWithNode<LO,GO,NO>(rows, comm);

    return Teuchos::rcp(new FortranOperator(callback, map, map));
  }

}
