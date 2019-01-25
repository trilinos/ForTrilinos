/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#ifndef FORTRILINOS_NOX_SOLVER_HPP
#define FORTRILINOS_NOX_SOLVER_HPP

#include <NOX.H>
#include <NOX_Thyra.H>
#include <NOX_Thyra_MatrixFreeJacobianOperator.hpp>
#include <NOX_MatrixFree_ModelEvaluatorDecorator.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

#include "model_evaluator.hpp"

namespace ForTrilinos {

  class NOXSolver {
    typedef double                                  SC;
    typedef int                                     LO;
    typedef long long                               GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;

    typedef Tpetra::MultiVector<SC,LO,GO,NO>        MultiVector;
    typedef ModelEvaluator<SC,LO,GO,NO>             ME;

  public:

    NOXSolver(const Teuchos::RCP<ME>& model);
    void setup(Teuchos::RCP<Teuchos::ParameterList>& plist);
    NOX::StatusTest::StatusType solve(Teuchos::RCP<MultiVector> initial_guess = Teuchos::null);
    Teuchos::RCP<MultiVector> get_solution() const;

  private:
    Teuchos::RCP<NOX::Solver::Generic> solver_;
    Teuchos::RCP<ME> model_;
  };
}

#endif
