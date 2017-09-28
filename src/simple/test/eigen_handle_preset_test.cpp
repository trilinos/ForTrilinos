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

    // Initialize MPI system
    Teuchos::GlobalMPISession mpiSession(&argc, &argv, NULL);
    RCP<const Teuchos::Comm<int> > comm = Teuchos::DefaultComm<int>::getComm();

    // Read in the parameter list
    Teuchos::RCP<ParameterList> paramList(new ParameterList());
    Teuchos::updateParametersFromXmlFileAndBroadcast("davidson.xml", Teuchos::Ptr<ParameterList>(paramList.get()), *comm);

    // For now, run without a preconditioner
    paramList->remove("Preconditioner Type", false/*throwIfExists*/);

    // Read in the matrices
    auto A = Tpetra::MatrixMarket::Reader<ForTrilinos::EigenHandle::Matrix>::readSparseFile("LHS_matrix.mat", comm);
    auto M = Tpetra::MatrixMarket::Reader<ForTrilinos::EigenHandle::Matrix>::readSparseFile("RHS_matrix.mat", comm);

    int n = A->getNodeNumRows();

    // The eigen solution
    std::vector<double> evectors(n);
    std::vector<double> evalues(1);

    // Step 1: initialize a handle
    ForTrilinos::EigenHandle si;
    si.init(comm);

    // Step 2: setup the problem
    si.setup_matrix(A);
    si.setup_matrix_rhs(M);

    // Step 3: setup the solver
    si.setup_solver(paramList);

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
