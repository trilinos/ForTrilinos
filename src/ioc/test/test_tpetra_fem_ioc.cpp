//@HEADER
// ************************************************************************
//
//            NOX: An Object-Oriented Nonlinear Solver Package
//                 Copyright (2002) Sandia Corporation
//
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger Pawlowski (rppawlo@sandia.gov) or
// Eric Phipps (etphipp@sandia.gov), Sandia National Laboratories.
// ************************************************************************
//  CVS Information
//  $Source$
//  $Author$
//  $Date$
//  $Revision$
// ************************************************************************
//@HEADER

// NOX Objects
#include "NOX.H"
#include "NOX_Thyra.H"

// Trilinos Objects
#include "Teuchos_Comm.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_XmlParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include "Tpetra_Core.hpp"
#include "Tpetra_Vector.hpp"

#include "BelosTypes.hpp"
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_Ifpack2PreconditionerFactory.hpp"

#include "NOX_Thyra_MatrixFreeJacobianOperator.hpp"
#include "NOX_MatrixFree_ModelEvaluatorDecorator.hpp"

#include "Tpetra_ModelEvaluator_1DFEM.hpp"

const Tpetra::global_size_t numGlobalElements = 100;


// Sets up and runs the nonlinear optimization in NOX
int main2(Teuchos::RCP<const Teuchos::Comm<int>>& comm,
          Teuchos::RCP<Teuchos::ParameterList>& plist)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  // Get default Tpetra template types
  using Scalar = Tpetra::Vector<>::scalar_type;
  using LO = Tpetra::Vector<>::local_ordinal_type;
  using GO = Tpetra::Vector<>::global_ordinal_type;
  using Node = Tpetra::Vector<>::node_type;

  RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

  Teuchos::TimeMonitor::zeroOutTimers();

  // Create the model evaluator object
  Scalar x00 = 0.0;
  Scalar x01 = 1.0;
  RCP<TpetraModelEvaluator1DFEM<Scalar,LO,GO,Node>> model =
    tpetraModelEvaluator1DFEM<Scalar,LO,GO,Node>(comm, numGlobalElements, x00, x01);

  Stratimikos::DefaultLinearSolverBuilder builder;

  {
    auto p = Teuchos::rcpFromRef(plist->sublist("Linear Solver Settings"));

    std::string prec = p->get("Preconditioner Type", "None");
    if (prec == "None") {
      // Do nothing
    } else if (prec == "Ifpack2") {
      using Base = Thyra::PreconditionerFactoryBase<Scalar>;
      using Impl = Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<Scalar,LO,GO,Node>>;
      builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "Preconditioner Type must be one of 'None', 'Ifpack2'")
    }
    builder.setParameterList(p);
  }

  RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar>> lows_factory =
    builder.createLinearSolveStrategy("");

  model->set_W_factory(lows_factory);

  // Create the initial guess
  RCP<Thyra::VectorBase<Scalar>> initial_guess =
    model->getNominalValues().get_x()->clone_v();
  Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<Scalar>::one());

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
      nox_group = rcp(new NOX::Thyra::Group(*initial_guess, model,
                                            model->create_W_op(), lows_factory,
                                            Teuchos::null, Teuchos::null));

    } else {
      auto p = jsettings.sublist("Jacobian Types");
      if (jtype == "Matrix Free Newton") {
        auto jfnk_params = Teuchos::rcpFromRef(p.sublist(jtype));
        Teuchos::ParameterList print_params;
        jfnk_op = rcp(new NOX::Thyra::MatrixFreeJacobianOperator<Scalar>(print_params));
        jfnk_op->setParameterList(jfnk_params);
        jfnk_params->print(*out);

        // Wrap the model evaluator in a JFNK Model Evaluator
        thyra_model = rcp(new NOX::MatrixFreeModelEvaluatorDecorator<Scalar>(model));

        // Create the Preconditioner operator
        if (jsettings.get<bool>("Use Prec", false))
          prec_op = thyra_model->create_W_prec();

        // Create the NOX::Thyra::Group
        nox_group = rcp(new NOX::Thyra::Group(*initial_guess, thyra_model, jfnk_op,
                                              lows_factory, prec_op, Teuchos::null));

      }
    }
  } else {
    // Default to Analytic Jacobian
    nox_group = rcp(new NOX::Thyra::Group(*initial_guess, model,
                                          model->create_W_op(), lows_factory,
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
  RCP<NOX::Solver::Generic> solver = NOX::Solver::buildSolver(nox_group, combo, nl_params);
  NOX::StatusTest::StatusType solvStatus = solver->solve();

  int returncode = (solvStatus == NOX::StatusTest::Converged) ? 0 : 1;

  Teuchos::TimeMonitor::summarize();

  return returncode;
}

int main(int argc, char* argv[])
{

  // GlobalMPISession calls MPI_Init() in its constructor, if
  // appropriate.
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);

  // Get a communicator corresponding to MPI_COMM_WORLD
  Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

  Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::createParameterList("Test Params");
  Teuchos::updateParametersFromXmlFileAndBroadcast("params.xml", params.ptr(), *comm);

  int errors = main2(comm, params);

  TEUCHOS_TEST_FOR_EXCEPTION(errors!=0,
   std::runtime_error,
   "One or more NOX evaluations failed!");

  if (comm->getRank() == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  return 0;

}
