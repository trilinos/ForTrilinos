/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#ifndef FORTRILINOS_SOLVER_HANDLE_HPP
#define FORTRILINOS_SOLVER_HANDLE_HPP

#include "ForTrilinos_config.h"

#include <Kokkos_DefaultNode.hpp>
#include "fortpetra/ForTrilinos_DefaultNodeType.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Thyra_LinearOpWithSolveBase.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <utility>

namespace ForTrilinos {

  class TrilinosSolver {
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
    typedef Thyra::LinearOpWithSolveBase<SC>        LOWS;

  public:

    // Constructors
    TrilinosSolver() : status_(NOT_INITIALIZED) { }

    TrilinosSolver(const TrilinosSolver&) = delete;
    void operator=(const TrilinosSolver&) = delete;

    // Initialize
    void init();
    void init(const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

    // Setup matrix
    void setup_matrix(const Teuchos::RCP<Matrix>& A);

    // Setup operator
    void setup_operator(const Teuchos::RCP<Operator>& A);

    // Setup solver based on the parameter list
    void setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    // Solve linear system given rhs
    void solve(const Teuchos::RCP<const MultiVector>& rhs, Teuchos::RCP<MultiVector>& lhs) const;

    // Free all data
    void finalize();

  private:
    Teuchos::RCP<const Teuchos::Comm<int>>  comm_;
    Teuchos::RCP<Operator>                  A_;
    Teuchos::RCP<ParameterList>             paramList_;
    Teuchos::RCP<LOWS>                      thyraInverseA_;

    enum Status {
      NOT_INITIALIZED = 0,
      INITIALIZED = 1,
      MATRIX_SETUP = 2,
      SOLVER_SETUP = 3,
    } status_;
  };

}
#endif // FORTRILINOS_SOLVER_HANDLE_HPP
