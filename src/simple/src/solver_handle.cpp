/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include "fortran_operator.hpp"
#include "solver_handle.hpp"
#include "handle_helpers.hpp"

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Teuchos_DefaultComm.hpp>

#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_IFPACK2
#  include <Teuchos_AbstractFactoryStd.hpp>
#  include <Thyra_Ifpack2PreconditionerFactory.hpp>
#endif

#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_MUELU
#  include <Stratimikos_MueLuHelpers.hpp>
#endif

#include <Thyra_TpetraLinearOp.hpp>

#include <stdexcept>

namespace ForTrilinos {

  void SolverHandle::init() {
    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);
    comm_ = Teuchos::DefaultComm<int>::getComm();
    status_ = INITIALIZED;
  }

  void SolverHandle::init(const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);
    comm_ = comm;
    status_ = INITIALIZED;
  }

  void SolverHandle::setup_matrix(std::pair<const GO*,size_t> rowInds, std::pair<const LO*,size_t> rowPtrs,
                                  std::pair<const GO*,size_t> colInds, std::pair<const SC*,size_t> values) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);
    auto A = HandleHelpers::setup_matrix_gen(comm_, rowInds, rowPtrs, colInds, values);
    setup_matrix(A);
  }

  void SolverHandle::setup_matrix(const Teuchos::RCP<Matrix> A) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);
    A_ = Teuchos::rcp_dynamic_cast<Operator>(A);
    status_ = MATRIX_SETUP;
  }

  void SolverHandle::setup_operator(std::pair<const GO*, size_t> rowInds, OperatorCallback callback) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);
    A_ = HandleHelpers::setup_operator_gen(comm_, rowInds, callback);
    status_ = MATRIX_SETUP;
  }

  void SolverHandle::setup_solver(const Teuchos::RCP<Teuchos::ParameterList> paramList) {
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
#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_IFPACK2
    {
      typedef Thyra::PreconditionerFactoryBase<SC>          Base;
      typedef Thyra::Ifpack2PreconditionerFactory<Matrix>   Impl;
      linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
    }
#endif
#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_MUELU
    Stratimikos::enableMueLu<LO,GO,NO>(linearSolverBuilder);
#endif

    linearSolverBuilder.setParameterList(paramList);

    // Build a new "solver factory" according to the previously specified parameter list
    auto solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);

    // Build a Thyra operator corresponding to A^{-1} computed using the Stratimikos solver
    thyraInverseA_ = Thyra::linearOpWithSolve(*solverFactory, A);

    status_ = SOLVER_SETUP;
  }

  void SolverHandle::solve(std::pair<const SC*, size_t> rhs, std::pair<SC*, size_t> lhs) const {
    auto map = A_->getDomainMap();
    auto size = lhs.second;

    TEUCHOS_ASSERT(size >= 0);
    TEUCHOS_ASSERT(lhs.second == size);
    TEUCHOS_ASSERT((rhs.first != NULL && lhs.second != NULL) || size == 0);
    TEUCHOS_ASSERT(map->getNodeNumElements() == size_t(size));

    // NOTE: This is a major simplification
    TEUCHOS_ASSERT(map->isSameAs(*(A_->getRangeMap())));

    auto X = Teuchos::rcp(new MultiVector(map, 1));
    auto B = Teuchos::rcp(new MultiVector(map, 1));

    // FIXME: data copying
    auto Xdata = X->getDataNonConst(0);
    auto Bdata = B->getDataNonConst(0);
    for (int i = 0; i < size; i++) {
      Xdata[i] = lhs.first[i];
      Bdata[i] = rhs.first[i];
    }

    solve(B, X);

    // FIXME: fix data copying
    for (int i = 0; i < size; i++)
      lhs.first[i] = Xdata[i];
  }

  void SolverHandle::solve(const Teuchos::RCP<const MultiVector> B, Teuchos::RCP<MultiVector> X) const {
    TEUCHOS_ASSERT(status_ == SOLVER_SETUP);

    TEUCHOS_ASSERT(X->getMap()->isSameAs(*A_->getDomainMap()));
    TEUCHOS_ASSERT(B->getMap()->isSameAs(*A_->getRangeMap()));

    auto map = A_->getDomainMap();

    auto thyraX = Thyra::createMultiVector(X);
    auto thyraB = Thyra::createConstMultiVector(B);
    // auto thyraX = Thyra::tpetraMultiVector<SC,LO,GO,NO>(Thyra::tpetraVectorSpace<SC,LO,GO,NO>(map), X);
    // auto thyraB = Thyra::tpetraMultiVector<SC,LO,GO,NO>(Thyra::tpetraVectorSpace<SC,LO,GO,NO>(map), B);

    auto status = Thyra::solve<SC>(*thyraInverseA_, Thyra::NOTRANS, *thyraB, thyraX.ptr());

    if (!map->getComm()->getRank())
      std::cout << status << std::endl;
  }

  void SolverHandle::finalize() {
    // No need to check the status_, we can finalize() at any moment.
    comm_          = Teuchos::null;
    A_             = Teuchos::null;
    paramList_     = Teuchos::null;
    thyraInverseA_ = Teuchos::null;

    status_ = NOT_INITIALIZED;
  }

}
