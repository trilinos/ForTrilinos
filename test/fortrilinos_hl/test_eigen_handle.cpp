/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include "fortrilinos_hl/eigen_handle.hpp"

#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <iostream>

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

    // Initialize MPI system
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // Read in the parameter list
    ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast("davidson.xml", Teuchos::Ptr<ParameterList>(&paramList), *comm);

    paramList.set("NumEV", 1);

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

    // The eigen solution
    std::vector<double> evalues(2);
    std::vector<int>    eindex(1);
    Teuchos::RCP<MultiVector> X = Teuchos::rcp(new MultiVector(rowMap, 1));

    // Step 1: initialize a handle
    ForTrilinos::TrilinosEigenSolver handle;
    handle.init(comm);

    // Step 2: setup the problem
    handle.setup_matrix(A);

    // Step 3: setup the solver
    handle.setup_solver(Teuchos::rcpFromRef(paramList));

    // Step 4: solve the system
    handle.solve(Teuchos::arrayViewFromVector(evalues), X, Teuchos::arrayViewFromVector(eindex));

    // TODO: Check the solution

    // Step 5: clean up
    handle.finalize();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
