/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include "eigen_handle.hpp"

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <AnasaziTpetraAdapter.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziFactory.hpp>

#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_IFPACK2
#include <Ifpack2_Factory.hpp>
#endif
#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_MUELU
#include <MueLu_CreateTpetraPreconditioner.hpp>
#endif

#include <Teuchos_TestForException.hpp>

#include <Thyra_TpetraLinearOp.hpp>

namespace ForTrilinos {

  void TrilinosEigenSolver::init() {
    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);
    comm_ = Teuchos::DefaultComm<int>::getComm();
    status_ = INITIALIZED;
  }

  void TrilinosEigenSolver::init(const Teuchos::RCP<const Teuchos::Comm<int>>& comm) {
    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);
    comm_ = comm;
    status_ = INITIALIZED;
  }

  void TrilinosEigenSolver::setup_matrix(const Teuchos::RCP<Matrix>& A) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);
    A_ = Teuchos::rcp_dynamic_cast<Operator>(A);
    status_ = MATRIX_SETUP;
  }
  void TrilinosEigenSolver::setup_matrix_rhs(const Teuchos::RCP<Matrix>& M) {
    TEUCHOS_ASSERT(status_ == INITIALIZED || status_ == MATRIX_SETUP);
    M_ = Teuchos::rcp_dynamic_cast<Operator>(M);
  }

  void TrilinosEigenSolver::setup_operator(OperatorCallback callback, const Teuchos::RCP<const Map>& domainMap, const Teuchos::RCP<const Map>& rangeMap) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);
    A_ = Teuchos::rcp(new FortranOperator(callback, domainMap, rangeMap));
    status_ = MATRIX_SETUP;
  }
  void TrilinosEigenSolver::setup_operator_rhs(OperatorCallback callback, const Teuchos::RCP<const Map>& domainMap, const Teuchos::RCP<const Map>& rangeMap) {
    TEUCHOS_ASSERT(status_ == INITIALIZED || status_ == MATRIX_SETUP);
    M_ = Teuchos::rcp(new FortranOperator(callback, domainMap, rangeMap));
  }

  void TrilinosEigenSolver::setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList) {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;

    TEUCHOS_ASSERT(status_ == MATRIX_SETUP);

    auto problem = Teuchos::rcp(new Anasazi::BasicEigenproblem<SC,MultiVector,Operator>());

    problem->setA(A_);
    problem->setM(M_);

    if (paramList->isParameter("Preconditioner Type")) {
      Teuchos::RCP<const Matrix> A = Teuchos::rcp_dynamic_cast<Matrix>(A_);
      TEUCHOS_TEST_FOR_EXCEPTION(A.is_null(), std::runtime_error,
        "Cannot set up a preconditioner for an implicit operator");

      auto prec_type = paramList->get<std::string>("Preconditioner Type");
      if (prec_type != "MueLu") {
#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_IFPACK2
        // Assume Ifpack2 type if not multigrid
        auto prec = Ifpack2::Factory::create(prec_type, A);
        if (paramList->isSublist(prec_type))
          prec->setParameters(paramList->sublist(prec_type));

        prec->initialize();
        prec->compute();

        problem->setPrec(prec);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(false, std::runtime_error,
            "Cannot create Ifpack2 preconditioner as Ifpack2 is not enabled");
#endif
      } else {
#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_MUELU
        Teuchos::ParameterList dummyList;
        auto prec = MueLu::CreateTpetraPreconditioner(A_,
            (paramList->isSublist(prec_type) ? paramList->sublist(prec_type) : dummyList));

        problem->setPrec(prec);
#else
        TEUCHOS_TEST_FOR_EXCEPTION(false, std::runtime_error,
            "Cannot create MueLu preconditioner as MueLu is not enabled");
#endif
      }
    }

    int numEV = paramList->get("NumEV", 1);
    problem->setNEV(numEV);

    // set random initial guess
    auto initVec = rcp(new MultiVector(A_->getDomainMap(), 1));
    problem->setInitVec(initVec);

    bool r = problem->setProblem();
    TEUCHOS_ASSERT(r);

    auto solverName = paramList->get<std::string>("Solver Type");
    solver_ = Anasazi::Factory::create(solverName, problem, paramList->sublist(solverName));

    status_ = SOLVER_SETUP;
  }

  void TrilinosEigenSolver::solve(std::pair<SC*, size_t> eigenValues, Teuchos::RCP<MultiVector>& eigenVectors) const {
    using Teuchos::RCP;
    using Teuchos::ArrayRCP;

    TEUCHOS_ASSERT(status_ == SOLVER_SETUP);

    Anasazi::ReturnType r = solver_->solve();
    TEUCHOS_ASSERT(r == 0);

    // Extract solution
    Anasazi::Eigensolution<SC,MultiVector> solution = solver_->getProblem().getSolution();

    std::vector<Anasazi::Value<SC>>& evalues = solution.Evals;
    for (int i = 0; i < std::min(eigenValues.second, evalues.size()); i++)
      eigenValues.first[i] = evalues[i].realpart; // FIXME: need to implement complex case

    eigenVectors = solution.Evecs;
  }

  void TrilinosEigenSolver::finalize() {
    // No need to check the status_, we can finalize() at any moment.
    comm_          = Teuchos::null;
    A_             = Teuchos::null;
    M_             = Teuchos::null;
    solver_        = Teuchos::null;
    paramList_     = Teuchos::null;

    status_ = NOT_INITIALIZED;
  }

}
