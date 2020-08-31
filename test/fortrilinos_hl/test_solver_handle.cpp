/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include "fortrilinos_hl/solver_handle.hpp"

#include <iostream>

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Teuchos_Array.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include "fortpetra/ForTrilinos_DefaultNodeType.hpp"

int main(int argc, char *argv[]) {
  bool success = false;
  bool verbose = true;

  try {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::ParameterList;
    using Teuchos::tuple;

    typedef double                                  SC;
    typedef int                                     LO;
    typedef long long                               GO;
    typedef ForTrilinos::DefaultNodeType            NO;

    using Map         = Tpetra::Map<LO,GO,NO>;
    using MultiVector = Tpetra::MultiVector<SC,LO,GO,NO>;
    using Matrix      = Tpetra::CrsMatrix<SC,LO,GO,NO>;

    using STS = Teuchos::ScalarTraits<SC>;

    // Initialize MPI system
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    int myRank   = comm->getRank();
    int numProcs = comm->getSize();

    // Read in the parameter list
    ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast("stratimikos.xml", Teuchos::Ptr<ParameterList>(&paramList), *comm);

    // Get tolerance from the parameter list
    const auto& solverList = paramList.sublist("Linear Solver Types").sublist(paramList.get<std::string>("Linear Solver Type"));
    const auto& belosList = solverList.sublist("Solver Types").sublist(solverList.get<std::string>("Solver Type"));
    auto tol = belosList.get<double>("Convergence Tolerance");

    // Set parameters
    const int numMyElements = 50;
    int numGlobalElements = numMyElements*comm->getSize();

    // Step 0: Construct tri-diagonal matrix
    RCP<Map>    rowMap = Teuchos::rcp(new Map(numGlobalElements, numMyElements, 0, comm));
    RCP<Matrix> A = Teuchos::rcp(new Matrix(rowMap, 3));
    for (LO lclRow = 0; lclRow < static_cast<LO>(numMyElements); ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);

      if (gblRow == 0) {
        // A(0, 0:1) = [2, -1]
        A->insertGlobalValues (gblRow,
                               tuple<GO> (gblRow, gblRow + 1),
                               tuple<SC> (2.0, -1.0));
      } else if (static_cast<Tpetra::global_size_t> (gblRow) == static_cast<Tpetra::global_size_t>(numGlobalElements - 1)) {
        // A(N-1, N-2:N-1) = [-1, 2]
        A->insertGlobalValues (gblRow,
                               tuple<GO> (gblRow - 1, gblRow),
                               tuple<SC> (-1.0, 2.0));
      } else {
        // A(i, i-1:i+1) = [-1, 2, -1]
        A->insertGlobalValues (gblRow,
                               tuple<GO> (gblRow - 1, gblRow, gblRow + 1),
                               tuple<SC> (-1.0, 2.0, -1.0));
      }
    }
    A->fillComplete();

    // The solution is iota
    RCP<MultiVector> rhs = Teuchos::rcp(new MultiVector(rowMap, 1));
    rhs->putScalar(0.0);
    auto data = rhs->getDataNonConst(0);
    if (myRank == 0)
      data[0] = -1.0;
    if (myRank == numProcs-1)
      data[numMyElements-1] = numGlobalElements;

    RCP<MultiVector> lhs = Teuchos::rcp(new MultiVector(rowMap, 1));
    lhs->randomize();

    // Calculate initial residual
    Teuchos::Array<typename STS::magnitudeType> norms(1);
    RCP<MultiVector> residual = Teuchos::rcp(new MultiVector(rowMap, 1, false));
    A->apply(*lhs, *residual, Teuchos::NO_TRANS, 1., 0.);
    residual->update(1., *rhs, -1.);
    residual->norm2(norms);

    auto r0 = norms[0];

    // =======================================
    // Use case 1: setup and solve the problem
    // =======================================

    // Step 1: initialize a handle
    ForTrilinos::TrilinosSolver handle;
    handle.init(comm);

    // Step 2: setup the problem
    handle.setup_matrix(A);

    // Step 3: setup the solver
    handle.setup_solver(Teuchos::rcpFromRef(paramList));

    // Step 4: solve the system
    handle.solve(rhs, lhs);

    // Check the solution
    A->apply(*lhs, *residual, Teuchos::NO_TRANS, 1., 0.);
    residual->update(1., *rhs, -1.);
    residual->norm2(norms);

    TEUCHOS_ASSERT(norms[0]/r0 < tol);

    // ========================================================
    // Use case 2: Reuse preconditioner with a different matrix
    // ========================================================
    A->resumeFill();
    for (LO lclRow = 0; lclRow < static_cast<LO>(numMyElements); ++lclRow) {
      const GO gblRow = rowMap->getGlobalElement (lclRow);

      A->replaceGlobalValues(gblRow, tuple<GO>(gblRow), tuple<SC>(3.0));
    }
    A->fillComplete();
    lhs->randomize();
    A->apply(*lhs, *residual, Teuchos::NO_TRANS, 1., 0.);
    residual->update(1., *rhs, -1.);
    residual->norm2(norms);
    r0 = norms[0];

    // Step 2a: setup the problem that has already been setup
    handle.setup_matrix(A);

    // Step 4a: solve the system with the original preconditioner
    handle.solve(rhs, lhs);

    // Check the solution
    A->apply(*lhs, *residual, Teuchos::NO_TRANS, 1., 0.);
    residual->update(1., *rhs, -1.);
    residual->norm2(norms);

    TEUCHOS_ASSERT(norms[0]/r0 < tol);


    // Step 5: clean up
    handle.finalize();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
