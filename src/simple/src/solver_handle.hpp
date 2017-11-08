#ifndef FORTRILINOS_SOLVER_HANDLE_HPP
#define FORTRILINOS_SOLVER_HANDLE_HPP

#include "ForTrilinosSimpleInterface_config.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Thyra_LinearOpWithSolveBase.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include "fortran_operator.hpp"

#include <utility>

namespace ForTrilinos {

  class SolverHandle {
  public:
    typedef double                                  SC;
    typedef int                                     LO;
    typedef int                                     GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
    typedef size_t                                  global_size_t;

    typedef Tpetra::Map<LO,GO,NO>                   Map;
    typedef Tpetra::Operator<SC,LO,GO,NO>           Operator;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO,false>    Matrix;
    typedef Tpetra::Vector<SC,LO,GO,NO,false>       Vector;
    typedef Tpetra::MultiVector<SC,LO,GO,NO,false>  MultiVector;
    typedef Teuchos::ParameterList                  ParameterList;
    typedef Thyra::LinearOpWithSolveBase<SC>        LOWS;

  public:

    // Constructors
    SolverHandle() : status_(NOT_INITIALIZED) { }

    SolverHandle(const SolverHandle&) = delete;
    void operator=(const SolverHandle&) = delete;

    // Initialize
    void init();
    void init(const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

    // Setup matrix
    void setup_matrix(std::pair<const int*,size_t> rowInds, std::pair<const int*,size_t> rowPtrs,
                      std::pair<const int*,size_t> colInds, std::pair<const double*,size_t> values);
    void setup_matrix(Teuchos::RCP<Matrix> A);

    // Setup operator
    void setup_operator(std::pair<const int*, size_t> rowInds, OperatorCallback callback);

    // Setup solver based on the parameter list
    void setup_solver(const Teuchos::RCP<Teuchos::ParameterList> paramList);

    // Solve linear system given rhs
    void solve(std::pair<const double*, size_t> rhs, std::pair<double*, size_t> lhs) const;
    void solve(const Teuchos::RCP<const MultiVector> rhs, Teuchos::RCP<MultiVector> lhs) const;

    // Free all data
    void finalize();

  private:
    Teuchos::RCP<const Teuchos::Comm<int>>  comm_;
    Teuchos::RCP<Operator>                  A_;
    Teuchos::RCP<ParameterList>             paramList_;
    Teuchos::RCP<LOWS>                      thyraInverseA_;

    enum Status {
      NOT_INITIALIZED,
      INITIALIZED,
      MATRIX_SETUP,
      SOLVER_SETUP,
    } status_;
  };

}
#endif // FORTRILINOS_SOLVER_HANDLE_HPP
