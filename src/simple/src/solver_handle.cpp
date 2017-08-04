#include "fortran_operator.hpp"
#include "solver_handle.hpp"
#include "handle_helpers.hpp"

#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <stdexcept>

#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_IFPACK2
#  include <Teuchos_AbstractFactoryStd.hpp>
#  include <Thyra_Ifpack2PreconditionerFactory.hpp>
#endif

#ifdef HAVE_FORTRILINOSSIMPLEINTERFACE_MUELU
#  include <Stratimikos_MueLuHelpers.hpp>
#endif

#include <Teuchos_DefaultMpiComm.hpp>
#include <Teuchos_DefaultSerialComm.hpp>

#include <Thyra_TpetraLinearOp.hpp>


namespace ForTrilinos {

  void SolverHandle::init() {
    using Teuchos::rcp;

    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);

    comm_ = rcp(new Teuchos::SerialComm<int>());

    status_ = INITIALIZED;
  }

  void SolverHandle::init(MPI_Comm comm) {
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

  void SolverHandle::setup_matrix(int numRows, const int* rowInds, const int* rowPtrs, int numNnz, const int* colInds, const double* values) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);

    A_ = HandleHelpers::setup_matrix_gen(comm_, numRows, rowInds, rowPtrs, numNnz, colInds, values);

    status_ = MATRIX_SETUP;
  }

  void SolverHandle::setup_operator(int numRows, const int* rowInds, OperatorCallback funcptr) {
    TEUCHOS_ASSERT(status_ == INITIALIZED);

    A_ = HandleHelpers::setup_operator_gen(comm_, numRows, rowInds, funcptr);

    status_ = MATRIX_SETUP;
  }

  void SolverHandle::setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList) {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_implicit_cast;

    TEUCHOS_ASSERT(status_ == MATRIX_SETUP);

    // FIXME: add validation check

    // Create a Thyra linear operator (A) using a Tpetra::CrsMatrix (A_)
    RCP<const Thyra::LinearOpBase<SC> > A = Thyra::tpetraLinearOp<SC,LO,GO,NO>(
          Thyra::tpetraVectorSpace<SC,LO,GO,NO>(A_->getDomainMap()),
          Thyra::tpetraVectorSpace<SC,LO,GO,NO>(A_->getRangeMap()),
          rcp_implicit_cast<Tpetra::Operator<SC,LO,GO,NO> >(A_));


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

  void SolverHandle::solve(int size, const double* rhs, double* lhs) const {
    using Teuchos::RCP;
    using Teuchos::ArrayRCP;

    TEUCHOS_ASSERT(status_ == SOLVER_SETUP);

    auto map = A_->getDomainMap();

    TEUCHOS_ASSERT(size >= 0);
    TEUCHOS_ASSERT((rhs != NULL && lhs != NULL) || size == 0);
    TEUCHOS_ASSERT(map->getNodeNumElements() == size_t(size));

    // NOTE: This is a major simplification
    TEUCHOS_ASSERT(map->isSameAs(*(A_->getRangeMap())));

    RCP<Vector> X = rcp(new Vector(map));
    RCP<Vector> B = rcp(new Vector(map));

    // FIXME: data copying
    auto Xdata = X->getDataNonConst();
    auto Bdata = B->getDataNonConst();
    for (int i = 0; i < size; i++) {
      Xdata[i] = lhs[i];
      Bdata[i] = rhs[i];
    }

    auto thyraX =
        Thyra::tpetraVector<SC,LO,GO,NO>(Thyra::tpetraVectorSpace<SC,LO,GO,NO>(map), X);
    auto thyraB =
        Thyra::tpetraVector<SC,LO,GO,NO>(Thyra::tpetraVectorSpace<SC,LO,GO,NO>(map), B);

    Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*thyraInverseA_, Thyra::NOTRANS, *thyraB, thyraX.ptr());

    if (!map->getComm()->getRank())
      std::cout << status << std::endl;

    // FIXME: fix data copying
    for (int i = 0; i < size; i++)
      lhs[i] = Xdata[i];
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
