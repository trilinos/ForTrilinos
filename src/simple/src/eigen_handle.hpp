/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#ifndef FORTRILINOS_EIGEN_HANDLE_HPP
#define FORTRILINOS_EIGEN_HANDLE_HPP

#include "ForTrilinosSimpleInterface_config.hpp"

#include <Teuchos_Comm.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <AnasaziSolverManager.hpp>

#include "fortran_operator.hpp"

#include <utility>

namespace ForTrilinos {

  class EigenHandle {
  public:
    typedef double                                  SC;
    typedef int                                     LO;
    typedef long long                               GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
    typedef size_t                                  global_size_t;

    typedef Teuchos::ParameterList                  ParameterList;
    typedef Tpetra::CrsMatrix<SC,LO,GO,NO>          Matrix;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>        MultiVector;
    typedef Tpetra::Operator<SC,LO,GO,NO>           Operator;
    typedef Tpetra::Vector<SC,LO,GO,NO>             Vector;
    typedef Anasazi::SolverManager<SC,MultiVector,Operator> SolverManager;

  public:

    // Constructors
    EigenHandle() : status_(NOT_INITIALIZED) { }

    EigenHandle(const EigenHandle&) = delete;
    void operator=(const EigenHandle&) = delete;

    // Initialize
    void init();
    void init(const Teuchos::RCP<const Teuchos::Comm<int>>& comm);

    // Setup matrices by construction
    void setup_matrix    (std::pair<const GO*,size_t> rowInds, std::pair<const LO*,size_t> rowPtrs,
                          std::pair<const GO*,size_t> colInds, std::pair<const SC*,size_t> values);
    void setup_matrix_rhs(std::pair<const GO*,size_t> rowInds, std::pair<const LO*,size_t> rowPtrs,
                          std::pair<const GO*,size_t> colInds, std::pair<const SC*,size_t> values);

    // Setup matrices by user provided data
    void setup_matrix(Teuchos::RCP<Matrix> A);
    void setup_matrix_rhs(Teuchos::RCP<Matrix> M);

    // Setup operators
    void setup_operator    (std::pair<const GO*, size_t> rowInds, OperatorCallback callback);
    void setup_operator_rhs(std::pair<const GO*, size_t> rowInds, OperatorCallback callback);

    // Setup solver based on the parameter list
    void setup_solver(const Teuchos::RCP<Teuchos::ParameterList>& paramList);

    // Solve linear system given rhs
    void solve(std::pair<SC*, size_t> eigenValues, std::pair<SC*, size_t> eigenVectors) const;

    // Free all data
    void finalize();

  private:

    Teuchos::RCP<const Teuchos::Comm<int>> comm_;
    Teuchos::RCP<Operator>           A_, M_;
    Teuchos::RCP<SolverManager>      solver_;
    Teuchos::RCP<ParameterList>      paramList_;

    enum Status {
      NOT_INITIALIZED,
      INITIALIZED,
      MATRIX_SETUP,
      SOLVER_SETUP,
    } status_;
  };

}
#endif // FORTRILINOS_EIGEN_HANDLE_HPP
