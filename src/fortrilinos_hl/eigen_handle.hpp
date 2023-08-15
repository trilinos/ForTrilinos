/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#ifndef FORTRILINOS_EIGEN_HANDLE_HPP
#define FORTRILINOS_EIGEN_HANDLE_HPP

#include "ForTrilinos_config.h"

#include <utility>

#include "fortpetra/ForTrilinos_DefaultNodeType.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <AnasaziSolverManager.hpp>


namespace ForTrilinos {

  class TrilinosEigenSolver {
  public:
    typedef double                                  SC;
    typedef int                                     LO;
    typedef long long                               GO;
    typedef ForTrilinos::DefaultNodeType            NO;
    typedef size_t                                  global_size_t;

    typedef Teuchos::ParameterList                  ParameterList;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO>          Matrix;
    typedef Tpetra::Map<LO,GO,NO>                   Map;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>        MultiVector;
    typedef Tpetra::Operator<SC,LO,GO,NO>           Operator;
    typedef Anasazi::SolverManager<SC,MultiVector,Operator> SolverManager;

  public:

    // Constructors
    TrilinosEigenSolver();
    explicit TrilinosEigenSolver(const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

    // Setup matrices by construction
    void setup_matrix(const Teuchos::RCP<Matrix>& A);
    void setup_matrix_rhs(const Teuchos::RCP<Matrix>& M);

    // Setup operators
    void setup_operator(const Teuchos::RCP<Operator>& A);
    void setup_operator_rhs(const Teuchos::RCP<Operator>& M);

    // Setup solver based on the parameter list
    void setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList);
    // Maximum number of eigenvalues to solve for
    int max_eigenvalues() const;

    // Solve eigen system given rhs
    int solve(Teuchos::ArrayView<SC> eigenValues,
              Teuchos::RCP<MultiVector>& eigenVectors,
              Teuchos::ArrayView<int> eigenIndex);
    // Whether the solver converged
    bool converged() const;
    // Number of iterations used to solve
    int num_iters() const;

  private:

    Teuchos::RCP<const Teuchos::Comm<int>> comm_;
    Teuchos::RCP<Operator>           A_, M_;
    Teuchos::RCP<SolverManager>      solver_;
    int                              numEigenvalues_;
    int                              numIters_;
    bool                             converged_;

    enum Status {
      INITIALIZED,
      MATRIX_SETUP,
      SOLVER_SETUP,
      SOLVED,
    } status_;
  };

}
#endif // FORTRILINOS_EIGEN_HANDLE_HPP
