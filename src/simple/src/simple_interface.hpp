#include <mpi.h>

#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Thyra_LinearOpWithSolveBase.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include "ForTrilinosSimpleInterface_config.hpp"

namespace ForTrilinos {

  class SimpleInterface {
  private:
    typedef double                                  SC;
    typedef int                                     LO;
    typedef int                                     GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
    typedef size_t                                  global_size_t;

    typedef Tpetra::Map<LO,GO,NO>                   Map;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO>          Matrix;
    typedef Tpetra::Vector<SC,LO,GO,NO>             Vector;
    typedef Teuchos::ParameterList                  ParameterList;
    typedef Thyra::LinearOpWithSolveBase<SC>        LOWS;

  private:
    // Constructors
    SimpleInterface() : status_(NOT_INITIALIZED) { }

  public:

    static SimpleInterface& getInstance() {
      static SimpleInterface instance;
      return instance;
    }

    SimpleInterface(const SimpleInterface&) = delete;
    void operator=(const SimpleInterface&)  = delete;

    // Initialize
    void init(MPI_Comm comm = MPI_COMM_WORLD);

    // Setup matrix
    void setup_matrix(int numRows, const int* rowInds, const int* rowPtrs, int numNnz, const int* colInds, const double* values);

    // Setup solver based on the parameter list
    void setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    // Solve linear system given rhs
    void solve(int size, const double* rhs, double* lhs) const;

    // Free all data
    void finalize();

  private:
    Teuchos::RCP<Teuchos::Comm<int>> comm_;
    Teuchos::RCP<Matrix>             A_;
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
