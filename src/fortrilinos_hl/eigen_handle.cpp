/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include "eigen_handle.hpp"
#include "fortrilinos_utilities.hpp"

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Teuchos_DefaultComm.hpp>

#include <AnasaziTpetraAdapter.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziFactory.hpp>

#if FORTRILINOS_USE_IFPACK2
#include <Ifpack2_Factory.hpp>
#endif
#if FORTRILINOS_USE_MUELU
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

  void TrilinosEigenSolver::setup_operator(const Teuchos::RCP<Operator>& A) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);
    A_ = A;
    status_ = MATRIX_SETUP;
  }
  void TrilinosEigenSolver::setup_operator_rhs(const Teuchos::RCP<Operator>& M) {
    TEUCHOS_ASSERT(status_ == INITIALIZED || status_ == MATRIX_SETUP);
    M_ = M;
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
#if FORTRILINOS_USE_IFPACK2
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
#if FORTRILINOS_USE_MUELU
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

    TEUCHOS_TEST_FOR_EXCEPTION(!paramList->isParameter("NumEV"), std::runtime_error,
        "Please specify the desired number of eigenvectors by providing \"NumEV\" parameter");
    numEigenvalues_ = paramList->get<int>("NumEV");
    problem->setNEV(numEigenvalues_);

    // set random initial guess
    auto initVec = rcp(new MultiVector(A_->getDomainMap(), 1));
    SetRandomSeed(*comm_);
    initVec->randomize();
    problem->setInitVec(initVec);

    bool r = problem->setProblem();
    TEUCHOS_ASSERT(r);

    TEUCHOS_TEST_FOR_EXCEPTION(!paramList->isParameter("Solver Type"), std::runtime_error,
        "Please specify the desired solver by providing \"Solver Type\" parameter");
    auto solverName = paramList->get<std::string>("Solver Type");
    solver_ = Anasazi::Factory::create(solverName, problem, paramList->sublist(solverName));

    status_ = SOLVER_SETUP;
  }

  int TrilinosEigenSolver::solve(Teuchos::ArrayView<SC> eigenValues,
                                 Teuchos::RCP<MultiVector>& eigenVectors,
                                 Teuchos::ArrayView<int> eigenIndex) const {
    using Teuchos::RCP;
    using Teuchos::ArrayRCP;

    TEUCHOS_ASSERT(status_ == SOLVER_SETUP);

    Anasazi::ReturnType r = solver_->solve();
    TEUCHOS_ASSERT(r == 0);

    // Extract solution
    Anasazi::Eigensolution<SC,MultiVector> solution = solver_->getProblem().getSolution();

    int eNum = solution.numVecs;
    TEUCHOS_ASSERT(eNum > 0);
    std::vector<Anasazi::Value<SC>>& eValues = solution.Evals;
    std::vector<int>&                eIndex  = solution.index;

    size_t numConverged = std::min(eNum, numEigenvalues_);
    TEUCHOS_TEST_FOR_EXCEPTION(2 * eigenValues.size() < numConverged, std::runtime_error,
      "Insufficient space to store eigenvalues. Please provide at least two times the desired number of eigenvalues.");
    TEUCHOS_TEST_FOR_EXCEPTION(eigenIndex.size() < numConverged, std::runtime_error,
      "Insufficient space to store index. Please provide at least two times the desired number of eigenvalues.");

    for (size_t i = 0; i < eValues.size(); i++) {
      eigenValues[2*i+0] = eValues[i].realpart;
      eigenValues[2*i+1] = eValues[i].imagpart;
    }

    for (size_t i = 0; i < eIndex.size(); i++)
      eigenIndex[i] = eIndex[i];

    eigenVectors = solution.Evecs;

    return eValues.size();
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
