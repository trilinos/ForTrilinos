/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include "solver_handle.hpp"

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Teuchos_DefaultComm.hpp>

#if FORTRILINOS_USE_IFPACK2
#  include <Teuchos_AbstractFactoryStd.hpp>
#  include <Thyra_Ifpack2PreconditionerFactory.hpp>
#endif

#if FORTRILINOS_USE_MUELU
#  include <Stratimikos_MueLuHelpers.hpp>
#endif

#include <Thyra_TpetraLinearOp.hpp>

#include <stdexcept>

namespace ForTrilinos {

  void TrilinosSolver::init() {
    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);
    comm_ = Teuchos::DefaultComm<int>::getComm();
    status_ = INITIALIZED;
  }

  void TrilinosSolver::init(const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);
    comm_ = comm;
    status_ = INITIALIZED;
  }

  void TrilinosSolver::setup_matrix(const Teuchos::RCP<Matrix>& A) {
    TEUCHOS_ASSERT(status_ >= INITIALIZED);
    A_ = Teuchos::rcp_dynamic_cast<Operator>(A);
    if (status_ < MATRIX_SETUP)
      status_ = MATRIX_SETUP;
  }

  void TrilinosSolver::setup_operator(const Teuchos::RCP<Operator>& A) {
    TEUCHOS_ASSERT(status_ >= INITIALIZED);
    A_ = A;
    if (status_ < MATRIX_SETUP)
      status_ = MATRIX_SETUP;
  }

  void TrilinosSolver::setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList) {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;

    TEUCHOS_ASSERT(status_ == MATRIX_SETUP);

    // FIXME: add validation check

    // Create a Thyra linear operator (A) using a Tpetra::CrsMatrix (A_)
    RCP<const Thyra::LinearOpBase<SC> > A = Thyra::tpetraLinearOp<SC,LO,GO,NO>(
          Thyra::tpetraVectorSpace<SC,LO,GO,NO>(A_->getDomainMap()),
          Thyra::tpetraVectorSpace<SC,LO,GO,NO>(A_->getRangeMap()),
          A_);


    // Build Stratimikos solver
    Stratimikos::DefaultLinearSolverBuilder linearSolverBuilder;
#if FORTRILINOS_USE_IFPACK2
    {
      typedef Thyra::PreconditionerFactoryBase<SC>          Base;
      typedef Thyra::Ifpack2PreconditionerFactory<Matrix>   Impl;
      linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
    }
#endif
#if FORTRILINOS_USE_MUELU
    Stratimikos::enableMueLu<LO,GO,NO>(linearSolverBuilder);
#endif

    linearSolverBuilder.setParameterList(paramList);

    // Build a new "solver factory" according to the previously specified parameter list
    auto solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);

    // Build a Thyra operator corresponding to A^{-1} computed using the Stratimikos solver
    thyraInverseA_ = Thyra::linearOpWithSolve(*solverFactory, A);

    status_ = SOLVER_SETUP;
  }

  void TrilinosSolver::solve(const Teuchos::RCP<const MultiVector>& B, Teuchos::RCP<MultiVector>& X) const {
    TEUCHOS_ASSERT(status_ == SOLVER_SETUP);

    TEUCHOS_ASSERT(X->getMap()->isSameAs(*A_->getDomainMap()));
    TEUCHOS_ASSERT(B->getMap()->isSameAs(*A_->getRangeMap()));

    auto map = A_->getDomainMap();

    auto thyraX = Thyra::createMultiVector(X);
    auto thyraB = Thyra::createConstMultiVector(B);

    auto status = Thyra::solve<SC>(*thyraInverseA_, Thyra::NOTRANS, *thyraB, thyraX.ptr());

    // FIXME
    if (!map->getComm()->getRank())
      std::cout << status << std::endl;
  }

  void TrilinosSolver::finalize() {
    // No need to check the status_, we can finalize() at any moment.
    comm_          = Teuchos::null;
    A_             = Teuchos::null;
    paramList_     = Teuchos::null;
    thyraInverseA_ = Teuchos::null;

    status_ = NOT_INITIALIZED;
  }

}
