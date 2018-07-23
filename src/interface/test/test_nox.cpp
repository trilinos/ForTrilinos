/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_AbstractFactoryStd.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>

#include "Tpetra_ModelEvaluator_1DFEM.hpp"
#include "nox_solver.hpp"


// Sets up and runs the nonlinear optimization in NOX
int main2(Teuchos::RCP<const Teuchos::Comm<int>>& comm,
          Teuchos::RCP<Teuchos::ParameterList>& plist)
{

  using Teuchos::RCP;
  using Teuchos::rcp;

  // Get default Tpetra template types
  using Scalar = Tpetra::CrsMatrix<>::scalar_type;
  using LO = Tpetra::CrsMatrix<>::local_ordinal_type;
  using GO = Tpetra::CrsMatrix<>::global_ordinal_type;
  using Node = Tpetra::CrsMatrix<>::node_type;

  Teuchos::TimeMonitor::zeroOutTimers();

  // Create the model evaluator object
  RCP<TpetraModelEvaluator1DFEM<Scalar,LO,GO,Node>> evaluator =
    rcp(new TpetraModelEvaluator1DFEM<Scalar,LO,GO,Node>(comm, 100, 0.0, 1.0));
  evaluator->setup(plist);

  ForTrilinos::NOXSolver<Scalar,LO,GO,Node> nox_solver(evaluator);
  nox_solver.setup(plist);
  NOX::StatusTest::StatusType solve_status = nox_solver.solve();

  int returncode = (solve_status == NOX::StatusTest::Converged) ? 0 : 1;

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
  Teuchos::updateParametersFromXmlFileAndBroadcast("nox_params.xml", params.ptr(), *comm);

  int errors = main2(comm, params);

  TEUCHOS_TEST_FOR_EXCEPTION(errors!=0,
   std::runtime_error,
   "One or more NOX evaluations failed!");

  if (comm->getRank() == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }
  return 0;
}
