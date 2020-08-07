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

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class NOXSolver {

  public:

    NOXSolver(const Teuchos::RCP<ModelEvaluator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& model);
    void setup(Teuchos::RCP<Teuchos::ParameterList>& plist);
    NOX::StatusTest::StatusType solve();

  private:
    Teuchos::RCP<NOX::Solver::Generic> solver_;
    Teuchos::RCP<ModelEvaluator<Scalar, LocalOrdinal, GlobalOrdinal, Node>> model_;
  };
}

#ifndef SWIG
#include "nox_solver_def.hpp"
#endif

#endif
