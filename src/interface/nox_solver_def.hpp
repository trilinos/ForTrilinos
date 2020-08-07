/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include <NOX.H>
#include <NOX_Thyra.H>
#include <NOX_Thyra_MatrixFreeJacobianOperator.hpp>
#include <NOX_MatrixFree_ModelEvaluatorDecorator.hpp>

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>

namespace ForTrilinos {

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  NOXSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  NOXSolver(const Teuchos::RCP<ModelEvaluator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>& model) :
    model_(model) {}

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void NOXSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  setup(Teuchos::RCP<Teuchos::ParameterList>& plist)
  {

    using Teuchos::RCP;
    using Teuchos::rcp;

    // Create the JFNK operator
    std::string jtype;
    RCP<NOX::Thyra::Group> nox_group;
    RCP<Thyra::PreconditionerBase<Scalar>> prec_op;
    RCP<Thyra::ModelEvaluator<Scalar>> thyra_model;
    RCP<NOX::Thyra::MatrixFreeJacobianOperator<Scalar>> jfnk_op;

    if (plist->isSublist("Jacobian Settings")) {
      auto jsettings = plist->sublist("Jacobian Settings");
      jtype = jsettings.get("Jacobian Type", "Analytic");
      if (jtype == "Analytic") {

        // Create the NOX::Thyra::Group
        nox_group = rcp(new NOX::Thyra::Group(*(model_->initial_guess), model_,
                                              model_->create_W_op(), model_->lows_factory,
                                              Teuchos::null, Teuchos::null));

      }
      else {
        auto p = jsettings.sublist("Jacobian Types");
        if (jtype == "Matrix Free Newton") {
          auto jfnk_params = Teuchos::rcpFromRef(p.sublist(jtype));
          Teuchos::ParameterList print_params;
          jfnk_op = rcp(new NOX::Thyra::MatrixFreeJacobianOperator<Scalar>(print_params));
          jfnk_op->setParameterList(jfnk_params);
          /*
          RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();
          jfnk_params->print(*out);
          */

          // Wrap the model evaluator in a JFNK Model Evaluator
          thyra_model = rcp(new NOX::MatrixFreeModelEvaluatorDecorator<Scalar>(model_));

          // Create the Preconditioner operator
          if (jsettings.get<bool>("Use Prec", false))
            prec_op = thyra_model->create_W_prec();

          // Create the NOX::Thyra::Group
          nox_group = rcp(new NOX::Thyra::Group(*(model_->initial_guess), thyra_model, jfnk_op,
                                                model_->lows_factory, prec_op, Teuchos::null));

        }
      }
    }
    else {
      // Default to Analytic Jacobian
      nox_group = rcp(new NOX::Thyra::Group(*(model_->initial_guess), model_,
                                            model_->create_W_op(), model_->lows_factory,
                                            Teuchos::null, Teuchos::null));
    }

    if (nox_group.is_null())
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unkown Jacobian type " << jtype);

    nox_group->computeF();

    // VERY IMPORTANT!!!  jfnk object needs base evaluation objects.
    // This creates a circular dependency, so use a weak pointer.
    if (!jfnk_op.is_null())
      jfnk_op->setBaseEvaluationToNOXGroup(nox_group.create_weak());

    // Create the NOX status tests and the solver
    // Create the convergence tests
    RCP<NOX::StatusTest::Combo> converged = rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::AND));
    RCP<NOX::StatusTest::NormF> absresid = rcp(new NOX::StatusTest::NormF(1.0e-8));
    RCP<NOX::StatusTest::NormWRMS> wrms = rcp(new NOX::StatusTest::NormWRMS(1.0e-2, 1.0e-8));
    converged->addStatusTest(absresid);
    converged->addStatusTest(wrms);

    RCP<NOX::StatusTest::Combo> combo = rcp(new NOX::StatusTest::Combo(NOX::StatusTest::Combo::OR));
    RCP<NOX::StatusTest::MaxIters> maxiters = rcp(new NOX::StatusTest::MaxIters(20));
    RCP<NOX::StatusTest::FiniteValue> fv = rcp(new NOX::StatusTest::FiniteValue);
    combo->addStatusTest(fv);
    combo->addStatusTest(converged);
    combo->addStatusTest(maxiters);

    // Create the solver
    auto nl_params = Teuchos::sublist(plist, "Nonlinear Solver Settings");
    solver_ = NOX::Solver::buildSolver(nox_group, combo, nl_params);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  NOX::StatusTest::StatusType
  NOXSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
  solve()
  {
    return solver_->solve();
  }
}
