/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <MatrixMarket_Tpetra.hpp>

#include <iostream>

#include "eigen_handle.hpp"

int main(int argc, char *argv[]) {
  bool success = false;
  bool verbose = true;

  try {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::ParameterList;

    using Matrix      = ForTrilinos::TrilinosEigenSolver::Matrix;
    using MultiVector = ForTrilinos::TrilinosEigenSolver::MultiVector;

    // Initialize MPI system
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // Read in the parameter list
    ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast("davidson.xml", Teuchos::Ptr<ParameterList>(&paramList), *comm);

    paramList.set("NumEV", 1);

    // For now, run without a preconditioner
    paramList.remove("Preconditioner Type", false/*throwIfExists*/);

    // Read in the matrices
    auto A = Tpetra::MatrixMarket::Reader<Matrix>::readSparseFile("LHS_matrix.mat", comm);
    auto M = Tpetra::MatrixMarket::Reader<Matrix>::readSparseFile("RHS_matrix.mat", comm);

    // The eigen solution
    std::vector<double> evalues(1);
    std::vector<int>    eindex(1);
    Teuchos::RCP<MultiVector> X = Teuchos::rcp(new MultiVector(A->getRowMap(), 1));

    // Step 1: initialize a handle
    ForTrilinos::TrilinosEigenSolver handle;
    handle.init(comm);

    // Step 2: setup the problem
    handle.setup_matrix(A);
    handle.setup_matrix_rhs(M);

    // Step 3: setup the solver
    handle.setup_solver(Teuchos::rcpFromRef(paramList));

    // Step 4: solve the system
    handle.solve(std::make_pair(evalues.data(), 1), X, std::make_pair(eindex.data(), 1));

    // TODO: Check the solution

    // Step 5: clean up
    handle.finalize();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
