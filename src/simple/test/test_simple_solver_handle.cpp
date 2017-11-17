#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

#include <iostream>

#include "solver_handle.hpp"

int main(int argc, char *argv[]) {
  bool success = false;
  bool verbose = true;

  try {
    using Teuchos::RCP;
    using Teuchos::rcp;
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::ParameterList;

    // Initialize MPI system
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    int myRank   = comm->getRank();
    int numProcs = comm->getSize();

    // Read in the parameter list
    ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast("stratimikos.xml", Teuchos::Ptr<ParameterList>(&paramList), *comm);

    // Set parameters
    const int n   = 50;
    int       nnz = 3*n;

    // Step 0: Construct tri-diagonal matrix, and rhs
    std::vector<int> rowPtrs(n+1);
    std::vector<long long> rowInds(n), colInds(nnz);
    std::vector<double> values(nnz);

    rowPtrs[0] = 0;
    int curPos = 0, offset = n * myRank;
    for (int i = 0; i < n; i++) {
      if (i || myRank) {
        colInds[curPos] = offset + i-1;
        values [curPos] = -1;
        curPos++;
      }
      colInds[curPos] = offset + i;
      values [curPos] = 2;
      curPos++;
      if (i != n-1 || myRank != numProcs-1) {
        colInds[curPos] = offset + i+1;
        values [curPos] = -1;
        curPos++;
      }
      rowPtrs[i+1] = curPos;

      rowInds[i] = offset + i;
    }
    nnz = curPos;

    colInds.resize(nnz);
    values. resize(nnz);

    // The solution lhs[i] = i
    std::vector<double> lhs(n), rhs(n);
    rhs[0]   = (myRank == 0          ?       -1.0 : 0.0);
    rhs[n-1] = (myRank == numProcs-1 ? offset + n : 0.0);

    // Step 1: initialize a handle
    ForTrilinos::SolverHandle si;
    si.init(comm);

    // Step 2: setup the problem
    si.setup_matrix(std::make_pair(rowInds.data(), rowInds.size()), std::make_pair(rowPtrs.data(), rowPtrs.size()),
                    std::make_pair(colInds.data(), colInds.size()), std::make_pair(values.data(), values.size()));

    // Step 3: setup the solver
    si.setup_solver(Teuchos::rcpFromRef(paramList));

    // Step 4: solve the system
    si.solve(std::make_pair(rhs.data(), rhs.size()), std::make_pair(lhs.data(), lhs.size()));

    // Check the solution
    double norm = 0.0;
    for (int i = 0; i < n; i++)
      norm += (lhs[i] - (offset+i))*(lhs[i] - (offset+i));
    norm = sqrt(norm);

    std::cout << "norm = " << norm << std::endl;

    // TODO: Get the tolerance out of the parameter list
    TEUCHOS_ASSERT(norm < 1e-6);

    // Step 5: clean up
    si.finalize();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
