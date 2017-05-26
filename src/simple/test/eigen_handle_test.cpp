#include <Teuchos_DefaultComm.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Teuchos_StandardCatchMacros.hpp>

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

    // Initialize MPI system
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();
    int myRank   = comm->getRank();
    int numProcs = comm->getSize();

#ifdef HAVE_MPI
    RCP<const Teuchos::MpiComm<int> > tmpic = rcp_dynamic_cast<const Teuchos::MpiComm<int> >(comm);
    TEUCHOS_TEST_FOR_EXCEPTION(tmpic.is_null(), std::runtime_error, "Cannot cast base Teuchos::Comm to Teuchos::MpiComm object.");
    RCP<const Teuchos::OpaqueWrapper<MPI_Comm> > rawMpiComm = tmpic->getRawMpiComm();
#endif

    // Read in the parameter list
    ParameterList paramList;
    Teuchos::updateParametersFromXmlFileAndBroadcast("davidson.xml", Teuchos::Ptr<ParameterList>(&paramList), *comm);

    // Set parameters
    const int n   = 50;
    int       nnz = 3*n;

    // Step 0: Construct tri-diagonal matrix
    std::vector<int> rowPtrs(n+1), rowInds(n), colInds(nnz);
    std::vector<double> values(nnz);

    rowPtrs[0] = 0;
    int curPos = 0, offset = n * myRank;
    for (int i = 0; i < n; i++) {
#if 1
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
#else
      colInds[curPos] = offset + i;
      values [curPos] = 1.0;
      curPos++;
#endif
    }
    nnz = curPos;

    colInds.resize(nnz);
    values. resize(nnz);

    // The eigen solution
    std::vector<double> evectors(n);
    std::vector<double> evalues(1);

    // Step 1: initialize a handle
    ForTrilinos::EigenHandle si;
#ifdef HAVE_MPI
    si.init(*rawMpiComm);
#else
    si.init();
#endif

    // Step 2: setup the problem
    si.setup_matrix(n, rowInds.data(), rowPtrs.data(), nnz, colInds.data(), values.data());

    // Step 3: setup the solver
    si.setup_solver(Teuchos::rcpFromRef(paramList));

    // Step 4: solve the system
    si.solve(1, evalues.data(), n, evectors.data());

    // TODO: Check the solution

    // Step 5: clean up
    si.finalize();

    success = true;
  }
  TEUCHOS_STANDARD_CATCH_STATEMENTS(verbose, std::cerr, success);

  return ( success ? EXIT_SUCCESS : EXIT_FAILURE );
}
