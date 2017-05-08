#include "fortran_operator.hpp"
#include "trilinos_handle.hpp"

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

  void TrilinosHandle::init() {
    using Teuchos::rcp;

    TEUCHOS_ASSERT(status_ == NOT_INITIALIZED);

    comm_ = rcp(new Teuchos::SerialComm<int>());

    status_ = INITIALIZED;
  }

  void TrilinosHandle::init(MPI_Comm comm) {
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

  void TrilinosHandle::setup_matrix(int numRows, const int* rowInds, const int* rowPtrs, int numNnz, const int* colInds, const double* values) {
    using Teuchos::ArrayRCP;
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ArrayView;

    TEUCHOS_ASSERT(status_ == INITIALIZED);

    TEUCHOS_ASSERT(numRows >= 0);
    TEUCHOS_ASSERT((rowInds != NULL && rowPtrs != NULL) || numRows == 0);
    TEUCHOS_ASSERT(numNnz >= 0);
    TEUCHOS_ASSERT((colInds != NULL && values != NULL) || numNnz == 0);

    ArrayView<const GO> rows(rowInds, numRows);
    RCP<const Map> rowMap = Tpetra::createNonContigMapWithNode<LO,GO,NO>(rows, comm_, {});

    RCP<Matrix> A = rcp(new Matrix(rowMap, 1));

    // TODO: Can we use setAllValues?
    for (int i = 0; i < numRows; i++) {
      ArrayView<const GO> cols(colInds + rowPtrs[i], rowPtrs[i+1] - rowPtrs[i]);
      ArrayView<const SC> vals(values  + rowPtrs[i], rowPtrs[i+1] - rowPtrs[i]);

      A->insertGlobalValues(rowInds[i], cols, vals);
    }

    A->fillComplete();

    A_ = A;

    status_ = MATRIX_SETUP;
  }

  void TrilinosHandle::setup_operator(int numRows, const int* rowInds, OperatorCallback funcptr) {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::ArrayView;

    TEUCHOS_ASSERT(status_ == INITIALIZED);

    TEUCHOS_ASSERT(numRows >= 0);
    TEUCHOS_ASSERT(rowInds != NULL || numRows == 0);

    // NOTE: we make a major assumption on the maps:
    //      rowMap == domainMap == rangeMap
    ArrayView<const GO> rows(rowInds, numRows);
    RCP<const Map> map = Tpetra::createNonContigMapWithNode<LO,GO,NO>(rows, comm_, {});

    A_ = rcp(new FortranOperator(funcptr, map, map));

    status_ = MATRIX_SETUP;
  }

  void TrilinosHandle::setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList) {
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
    RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > solverFactory =
      Thyra::createLinearSolveStrategy(linearSolverBuilder);

    // Build a Thyra operator corresponding to A^{-1} computed using the Stratimikos solver
    thyraInverseA_ = Thyra::linearOpWithSolve(*solverFactory, A);

    status_ = SOLVER_SETUP;
  }

  void TrilinosHandle::solve(int size, const double* rhs, double* lhs) const {
    using Teuchos::RCP;
    using Teuchos::ArrayRCP;

    TEUCHOS_ASSERT(status_ == SOLVER_SETUP);

    RCP<const Map> map = A_->getDomainMap();

    TEUCHOS_ASSERT(size >= 0);
    TEUCHOS_ASSERT((rhs != NULL && lhs != NULL) || size == 0);
    TEUCHOS_ASSERT(map->getNodeNumElements() == size_t(size));

    // NOTE: This is a major simplification
    TEUCHOS_ASSERT(map->isSameAs(*(A_->getRangeMap())));

    RCP<Vector> X = rcp(new Vector(map));
    RCP<Vector> B = rcp(new Vector(map));

    // FIXME: data copying
    ArrayRCP<SC> Xdata = X->getDataNonConst();
    ArrayRCP<SC> Bdata = B->getDataNonConst();
    for (int i = 0; i < size; i++) {
      Xdata[i] = lhs[i];
      Bdata[i] = rhs[i];
    }

    RCP<      Thyra::VectorBase<SC> >thyraX = Thyra::tpetraVector<SC,LO,GO,NO>(Thyra::tpetraVectorSpace<SC,LO,GO,NO>(map), X);
    RCP<const Thyra::VectorBase<SC> >thyraB = Thyra::tpetraVector<SC,LO,GO,NO>(Thyra::tpetraVectorSpace<SC,LO,GO,NO>(map), B);

    Thyra::SolveStatus<SC> status = Thyra::solve<SC>(*thyraInverseA_, Thyra::NOTRANS, *thyraB, thyraX.ptr());

    if (!map->getComm()->getRank())
      std::cout << status << std::endl;

    // FIXME: fix data copying
    for (int i = 0; i < size; i++)
      lhs[i] = Xdata[i];
  }

  void TrilinosHandle::finalize() {
    // No need to check the status_, we can finalize() at any moment.
    comm_          = Teuchos::null;
    A_             = Teuchos::null;
    paramList_     = Teuchos::null;
    thyraInverseA_ = Teuchos::null;

    status_ = NOT_INITIALIZED;
  }

}
