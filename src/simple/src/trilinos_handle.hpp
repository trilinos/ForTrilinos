#ifndef FORTRILINOS_TRILINOS_HANDLE_HPP
#define FORTRILINOS_TRILINOS_HANDLE_HPP

#include "ForTrilinosSimpleInterface_config.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#else
typedef int MPI_Comm;
#endif

#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Thyra_LinearOpWithSolveBase.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include "fortran_operator.hpp"

namespace ForTrilinos {

  class TrilinosHandle {
  private:
    typedef double                                  SC;
    typedef int                                     LO;
    typedef int                                     GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
    typedef size_t                                  global_size_t;

    typedef Tpetra::Map<LO,GO,NO>                   Map;
    typedef Tpetra::Operator<SC,LO,GO,NO>           Operator;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO>          Matrix;
    typedef Tpetra::Vector<SC,LO,GO,NO>             Vector;
    typedef Teuchos::ParameterList                  ParameterList;
    typedef Thyra::LinearOpWithSolveBase<SC>        LOWS;
  public:
    typedef void (*OperatorCallback)(int n, const double* x, double* y);

  public:

    // Constructors
    TrilinosHandle() : status_(NOT_INITIALIZED) { }

    TrilinosHandle(const TrilinosHandle&) = delete;
    void operator=(const TrilinosHandle&) = delete;

    // Initialize
    void init();
    void init(MPI_Comm comm);

    // Setup matrix
    void setup_matrix(int numRows, const int* rowInds, const int* rowPtrs, int numNnz, const int* colInds, const double* values);

    // Setup operator
    void setup_operator(int numRows, const int* rowInds, OperatorCallback callback);

    // Setup solver based on the parameter list
    void setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    // Solve linear system given rhs
    void solve(int size, const double* rhs, double* lhs) const;

    // Free all data
    void finalize();

  private:
    Teuchos::RCP<Teuchos::Comm<int>> comm_;
    Teuchos::RCP<Operator>           A_;
    Teuchos::RCP<ParameterList>      paramList_;
    Teuchos::RCP<LOWS>               thyraInverseA_;

    enum Status {
      NOT_INITIALIZED,
      INITIALIZED,
      MATRIX_SETUP,
      SOLVER_SETUP,
    } status_;
  };

}
#endif // FORTRILINOS_TRILINOS_HANDLE_HPP
