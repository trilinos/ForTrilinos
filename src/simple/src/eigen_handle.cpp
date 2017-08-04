#include "eigen_handle.hpp"
#include "handle_helpers.hpp"

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>

#include <AnasaziTpetraAdapter.hpp>
#include <AnasaziBasicEigenproblem.hpp>
#include <AnasaziFactory.hpp>

#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_IFPACK2
#include <Ifpack2_Factory.hpp>
#endif
#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_MUELU
#include <MueLu_CreateTpetraPreconditioner.hpp>
#endif

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>
#include <Teuchos_TestForException.hpp>

#include <Thyra_TpetraLinearOp.hpp>

namespace ForTrilinos {

  void EigenHandle::init() {
    using Teuchos::rcp;

    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);

    comm_ = rcp(new Teuchos::SerialComm<int>());

    status_ = INITIALIZED;
  }

  void EigenHandle::init(MPI_Comm comm) {
#ifdef HAVE_MPI
    using Teuchos::rcp;

    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);

    {
      // Test if the communicator is valid
      int rank;
      int r = MPI_Comm_rank(comm, &rank);
      TEUCHOS_ASSERT(r == 0);
    }

    comm_ = rcp(new Teuchos::MpiComm<int>(comm));

    status_ = INITIALIZED;
#else
    throw std::runtime_error("MPI is not enabled");
#endif
  }

  void EigenHandle::setup_matrix(int numRows, const int* rowInds, const int* rowPtrs, int numNnz, const int* colInds, const double* values) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);

    A_ = HandleHelpers::setup_matrix_gen(comm_, numRows, rowInds, rowPtrs, numNnz, colInds, values);

    status_ = MATRIX_SETUP;
  }
  void EigenHandle::setup_matrix_rhs(int numRows, const int* rowInds, const int* rowPtrs, int numNnz, const int* colInds, const double* values) {
    TEUCHOS_ASSERT(status_ == INITIALIZED || status_ == MATRIX_SETUP);

    M_ = HandleHelpers::setup_matrix_gen(comm_, numRows, rowInds, rowPtrs, numNnz, colInds, values);
  }
  void EigenHandle::setup_operator(int numRows, const int* rowInds, void (*funcptr)(int n, const double* x, double* y)) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);

    A_ = HandleHelpers::setup_operator_gen(comm_, numRows, rowInds, funcptr);

    status_ = MATRIX_SETUP;
  }
  void EigenHandle::setup_operator_rhs(int numRows, const int* rowInds, void (*funcptr)(int n, const double* x, double* y)) {
    TEUCHOS_ASSERT(status_ == INITIALIZED || status_ == MATRIX_SETUP);

    M_ = HandleHelpers::setup_operator_gen(comm_, numRows, rowInds, funcptr);
  }

  void EigenHandle::setup_matrix(Teuchos::RCP<Matrix> A) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);

    A_ = A;

    status_ = MATRIX_SETUP;
  }
  void EigenHandle::setup_matrix_rhs(Teuchos::RCP<Matrix> M) {
    TEUCHOS_ASSERT(status_ == INITIALIZED || status_ == MATRIX_SETUP);

    M_ = M;
  }

  void EigenHandle::setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList) {
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
    TEUCHOS_ASSERT(numEV == 1);  // FIXME later
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

  void EigenHandle::solve(int numEigs, double* eigenValues, int size, double* eigenVectors) const {
    using Teuchos::RCP;
    using Teuchos::ArrayRCP;

    TEUCHOS_ASSERT(status_ == SOLVER_SETUP);

    Anasazi::ReturnType r = solver_->solve();
    TEUCHOS_ASSERT(r == 0);

    // Extract solution
    Anasazi::Eigensolution<SC,MultiVector> solution = solver_->getProblem().getSolution();

    std::vector<Anasazi::Value<SC>>& evalues = solution.Evals;
    RCP<MultiVector>& evectors = solution.Evecs;

    eigenValues[0] = evalues[0].realpart;

    memcpy(eigenVectors, evectors->getData(0).getRawPtr(), size*sizeof(double));
  }

  void EigenHandle::finalize() {
    // No need to check the status_, we can finalize() at any moment.
    comm_          = Teuchos::null;
    A_             = Teuchos::null;
    M_             = Teuchos::null;
    solver_        = Teuchos::null;
    paramList_     = Teuchos::null;

    status_ = NOT_INITIALIZED;
  }

}
