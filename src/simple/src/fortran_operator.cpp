#include "fortran_operator.hpp"

namespace ForTrilinos {

  void FortranOperator::apply(const MultiVector& X, MultiVector& Y,
                              Teuchos::ETransp mode, SC alpha, SC beta) const {
    using Teuchos::RCP;
    using Teuchos::ArrayRCP;
    typedef Teuchos::ScalarTraits<SC> STS;

    TEUCHOS_ASSERT(alpha == STS::one());
    TEUCHOS_ASSERT(beta  == STS::zero());

    TEUCHOS_ASSERT(X.getNumVectors() == 1);
    TEUCHOS_ASSERT(Y.getNumVectors() == 1);
    TEUCHOS_ASSERT(Y.getLocalLength() == X.getLocalLength());
    TEUCHOS_ASSERT(domainMap_->getNodeNumElements() == X.getLocalLength());
    TEUCHOS_ASSERT(rangeMap_ ->getNodeNumElements() == Y.getLocalLength());

    ArrayRCP<const SC> Xdata = X.getData(0);
    ArrayRCP<SC>       Ydata = Y.getDataNonConst(0);

    int n = Xdata.size();
    funcptr_(n, Xdata.get(), Ydata.get());
  }
}
