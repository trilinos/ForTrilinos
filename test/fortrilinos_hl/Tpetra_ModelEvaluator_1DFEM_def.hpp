/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#ifndef FORTRILINOS_TPETRA_MODELEVALUATOR_1DFEM_DEF_HPP
#define FORTRILINOS_TPETRA_MODELEVALUATOR_1DFEM_DEF_HPP

// Kokkos support
#include <Kokkos_Core.hpp>

// Nonmember constuctors

// Constructor

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
TpetraModelEvaluator1DFEM<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
TpetraModelEvaluator1DFEM(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                          const Tpetra::global_size_t num_global_elems,
                          const Scalar z_min,
                          const Scalar z_max) :
  comm_(comm)
{

  TEUCHOS_ASSERT(nonnull(comm_));

  const Tpetra::global_size_t num_nodes = num_global_elems + 1;

  // owned space
  GO idx_base = 0;
  x_owned_map_ = Teuchos::rcp(new const Map(num_nodes, idx_base, comm_));

  // ghosted space
  if (comm_->getSize() == 1) {
    x_ghosted_map_ = x_owned_map_;
  }
  else {
    GO min_overlap_GID;
    size_t num_overlap_nodes = x_owned_map_->getNodeNumElements() + 2;
    if ((comm_->getRank() == 0) || (comm_->getRank() == (comm_->getSize() - 1)))
      --num_overlap_nodes;

    if (comm_->getRank() == 0)
      min_overlap_GID = x_owned_map_->getMinGlobalIndex();
    else
      min_overlap_GID = x_owned_map_->getMinGlobalIndex() - 1;

    Teuchos::Array<GO> node_gids(num_overlap_nodes);
    GO gid = min_overlap_GID;
    for (auto node_gid=node_gids.begin(); node_gid!=node_gids.end(); ++node_gid) {
      *node_gid = gid;
      ++gid;
    }

    const Tpetra::global_size_t invalid = Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid();
    x_ghosted_map_ = Teuchos::rcp(new Map(invalid, node_gids, idx_base, comm_));
  }

  importer_ = Teuchos::rcp(new Tpetra::Import<LO,GO,NO>(x_owned_map_, x_ghosted_map_));

  // residual space
  f_owned_map_ = x_owned_map_;

  // Initialize the graph for W CrsMatrix object
  graph_ = create_graph(x_owned_map_, x_ghosted_map_);

  // Create the nodal coordinates
  node_coords_ = create_mesh(x_owned_map_, z_min, z_max, num_global_elems);

  // Timers
  resid_timer_ = Teuchos::TimeMonitor::getNewCounter("Model Evaluator: Residual Evaluation");
  jac_timer_ = Teuchos::TimeMonitor::getNewCounter("Model Evaluator: Jacobian Evaluation");

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinal, GlobalOrdinal, Node>>
TpetraModelEvaluator1DFEM<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
create_graph(const Teuchos::RCP<const Map>& owned_map,
             const Teuchos::RCP<const Map>& ghosted_map)
{

  // Create the shell for the graph
  Teuchos::RCP<Graph> graph =
    Teuchos::rcp(new Graph(owned_map, ghosted_map, 5));

  // Declare required variables
  GO row, column;
  size_t num_my_overlap_nodes = ghosted_map->getNodeNumElements();

  // Loop Over # of Finite Elements on Processor
  for (LO ne=0; ne<static_cast<LO>(num_my_overlap_nodes-1); ne++) {

    // Loop over nodes in element
    for (LO i=0; i<2; i++) {

      row = ghosted_map->getGlobalElement(ne+i);

      // Loop over Trial Functions
      for(LO j=0; j<2; j++) {

        // If this row is owned by current processor, add the index
        if (owned_map->isNodeGlobalElement(row)) {
          column = ghosted_map->getGlobalElement(ne+j);
          graph->insertGlobalIndices(row, 1, &column);
        }

      }
    }
  }
  graph->fillComplete();
  return graph;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Tpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
TpetraModelEvaluator1DFEM<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
create_mesh(const Teuchos::RCP<const Map>& owned_map,
            const SC z_min, const SC z_max,
            const Tpetra::global_size_t num_elems)
{
  size_t num_local_nodes = owned_map->getNodeNumElements();
  GO min_GID = owned_map->getMinGlobalIndex();
  SC dz = (z_max - z_min)/static_cast<SC>(num_elems);

  Teuchos::RCP<MultiVector> coords = Teuchos::rcp(new MultiVector(owned_map, 1));
  for (LO node=0; node<static_cast<LO>(num_local_nodes); node++) {
    coords->replaceLocalValue(node, 0, z_min + dz * static_cast<SC>(min_GID+node));
  }
  return coords;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<Tpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
TpetraModelEvaluator1DFEM<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
create_operator() const
{
  Teuchos::RCP<Operator> op = Teuchos::rcp(new Matrix(graph_));
  return op;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraModelEvaluator1DFEM<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
update_solution_vector(const Teuchos::RCP<const MultiVector>& xp) const
{
  typedef Kokkos::HostSpace host_space;

  // Create ghosted objects
  if (is_null(u_ptr_))
    u_ptr_ = Teuchos::rcp(new MultiVector(x_ghosted_map_, 1));
  u_ptr_->template modify<host_space>();
  u_ptr_->doImport(*xp, *importer_, Tpetra::REPLACE);

  if (is_null(x_ptr_)) {
    x_ptr_ = Teuchos::rcp(new MultiVector(x_ghosted_map_, 1));
    x_ptr_->template modify<host_space>();
    x_ptr_->doImport(*node_coords_, *importer_, Tpetra::INSERT);
  }

  x_ptr_->template sync<host_space>();
  u_ptr_->template sync<host_space>();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraModelEvaluator1DFEM<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
evaluate_residual(const Teuchos::RCP<const MultiVector>& xp,
                  Teuchos::RCP<MultiVector>& f) const
{

  // Update the solution variable
  update_solution_vector(xp);

  typedef Kokkos::HostSpace host_space;

  f->putScalar(0.0);
  f->template sync<host_space>();
  f->template modify<host_space>();

  Teuchos::TimeMonitor timer(*resid_timer_);
  const LO invalid = Tpetra::Details::OrdinalTraits<LO>::invalid();
  int my_rank = comm_->getRank();
  auto num_my_elems = x_ghosted_map_->getNodeNumElements()-1;

  // Loop Over # of Finite Elements on Processor
  auto x = x_ptr_->template getLocalView<Kokkos::HostSpace>();
  auto u = u_ptr_->template getLocalView<Kokkos::HostSpace>();

  for (LO ne=0; ne<static_cast<LO>(num_my_elems); ne++) {

    // Get the solution and coordinates at the nodes
    SC xx[2];
    xx[0] = x(ne, 0);
    xx[1] = x(ne+1, 0);

    SC uu[2];
    uu[0] = u(ne, 0);
    uu[1] = u(ne+1, 0);

    Linear2NodeFEBasis<SC, LO> basis;

    // Loop Over Gauss Points
    for(LO gp=0; gp < 2; gp++) {

      // Calculate the basis function at the gauss point
      basis.compute_basis(gp, xx, uu);

      // Loop over nodes in Element
      for (LO i=0; i< 2; i++) {
        LO local_row =
          x_owned_map_->getLocalElement(x_ghosted_map_->getGlobalElement(ne+i));
        if (local_row != invalid) {
          SC value = basis.wt * basis.dz *
            (basis.uu * basis.uu * basis.phi[i] + (basis.duu * basis.dphide[i])/(basis.dz * basis.dz));
          f->sumIntoLocalValue(local_row, 0, value);
        }
      }
    }

    // Correct for Dirichlet BCs
    if ((my_rank == 0) && (ne == 0)) {
      SC value = u(0, 0) - 1.0;
      f->replaceLocalValue(0, 0, value);
    }

  }

}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraModelEvaluator1DFEM<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
evaluate_jacobian(const Teuchos::RCP<const MultiVector>& xp,
                  Teuchos::RCP<Operator>& op) const
{

  // Update the solution variable
  update_solution_vector(xp);

  Teuchos::RCP<Matrix> J = Teuchos::rcp_dynamic_cast<Matrix>(op);
  TEUCHOS_ASSERT(nonnull(J));
  J->resumeFill();
  J->setAllToScalar(0.0);

  Teuchos::TimeMonitor timer(*jac_timer_);

  const LO invalid = Tpetra::Details::OrdinalTraits<LO>::invalid();
  int my_rank = comm_->getRank();
  auto num_my_elems = x_ghosted_map_->getNodeNumElements()-1;

  // Loop Over # of Finite Elements on Processor
  auto x = x_ptr_->template getLocalView<Kokkos::HostSpace>();
  auto u = u_ptr_->template getLocalView<Kokkos::HostSpace>();

  for (LO ne=0; ne<static_cast<LO>(num_my_elems); ne++) {

    // Get the solution and coordinates at the nodes
    SC xx[2];
    xx[0] = x(ne, 0);
    xx[1] = x(ne+1, 0);

    SC uu[2];
    uu[0] = u(ne, 0);
    uu[1] = u(ne+1, 0);

    Linear2NodeFEBasis<SC, LO> basis;

    // Loop Over Gauss Points
    for(LO gp=0; gp < 2; gp++) {

      // Calculate the basis function at the gauss point
      basis.compute_basis(gp, xx, uu);

      // Loop over NOs in Element
      for (LO i=0; i< 2; i++) {
        LO local_row =
          x_owned_map_->getLocalElement(x_ghosted_map_->getGlobalElement(ne+i));
        if (local_row != invalid) {
          // Loop over trial functions
          for (LO j=0; j<2; ++j) {
            const LO local_col = ne + j;
            SC value = basis.wt * basis.dz
              * ((basis.dphide[j]*basis.dphide[i])/(basis.dz*basis.dz)
              + 2.0*basis.uu*basis.phi[j]*basis.phi[i]);
            J->sumIntoLocalValues(local_row, 1, &value, &local_col);
          }
        }
      }
    }

    // Correct for Dirichlet BCs
    if ((my_rank == 0) && (ne == 0)) {
      LO row = 0;
      LO column = 0;
      SC value = 1.0;
      J->replaceLocalValues(row, 1, &value, &column);
      column = 1;
      value = 0.0;
      J->replaceLocalValues(row, 1, &value, &column);
    }
  }

  J->fillComplete();
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TpetraModelEvaluator1DFEM<Scalar, LocalOrdinal, GlobalOrdinal, Node>::
evaluate_preconditioner(const Teuchos::RCP<const MultiVector>& xp,
                        Teuchos::RCP<Operator>& op) const
{

  // Update the solution variable
  update_solution_vector(xp);

  Teuchos::RCP<Matrix> M_inv = Teuchos::rcp_dynamic_cast<Matrix>(op);
  TEUCHOS_ASSERT(nonnull(M_inv));
  M_inv->resumeFill();

  const LO invalid = Tpetra::Details::OrdinalTraits<LO>::invalid();
  int my_rank = comm_->getRank();
  auto num_my_elems = x_ghosted_map_->getNodeNumElements()-1;
  auto row_map = M_inv->getRowMap();
  auto col_map = M_inv->getColMap();

  if (is_null(J_diagonal_))
    J_diagonal_ = Teuchos::rcp(new Vector(x_owned_map_));
  M_inv->setAllToScalar(0.0);
  J_diagonal_->putScalar(0.0);

  Teuchos::TimeMonitor timer(*jac_timer_);

  // Loop Over # of Finite Elements on Processor
  auto x = x_ptr_->template getLocalView<Kokkos::HostSpace>();
  auto u = u_ptr_->template getLocalView<Kokkos::HostSpace>();

  for (LO ne=0; ne<static_cast<LO>(num_my_elems); ne++) {

    // Get the solution and coordinates at the nodes
    SC xx[2];
    xx[0] = x(ne, 0);
    xx[1] = x(ne+1, 0);

    SC uu[2];
    uu[0] = u(ne, 0);
    uu[1] = u(ne+1, 0);

    Linear2NodeFEBasis<SC, LO> basis;

    // Loop Over Gauss Points
    for(LO gp=0; gp < 2; gp++) {

      // Calculate the basis function at the gauss point
      basis.compute_basis(gp, xx, uu);

      // Loop over NOs in Element
      for (LO i=0; i< 2; i++) {
        LO local_row =
          x_owned_map_->getLocalElement(x_ghosted_map_->getGlobalElement(ne+i));
        if (local_row != invalid) {
          // Loop over trial functions
          for (LO j=0; j<2; ++j) {
            const LO local_col = ne + j;
            if (row_map->getGlobalElement(local_row) == col_map->getGlobalElement(local_col)) {
              SC value = basis.wt * basis.dz
                * ((basis.dphide[j]*basis.dphide[i])/(basis.dz*basis.dz)
                + 2.0*basis.uu*basis.phi[j]*basis.phi[i]);
              M_inv->sumIntoLocalValues(local_row, 1, &value, &local_col);
            }
          }
        }
      }
    }

    // Correct for Dirichlet BCs
    if ((my_rank == 0) && (ne == 0)) {
      LO row = 0;
      LO column = 0;
      SC value = 1.0;
      M_inv->replaceLocalValues(row, 1, &value, &column);
    }
  }

  // Invert the Jacobian diagonal for the preconditioner
  // For some reason the matrix must be fill complete before calling rightScale
  M_inv->fillComplete();
  Vector& diag = *J_diagonal_;
  M_inv->getLocalDiagCopy(diag);
  diag.reciprocal(diag);
  M_inv->rightScale(diag);
  M_inv->rightScale(diag);
}

// Finite Element Basis Object
template <class Scalar, class LocalOrdinal>
Linear2NodeFEBasis<Scalar, LocalOrdinal>::Linear2NodeFEBasis():
  uu(0.0),
  zz(0.0),
  duu(0.0),
  eta(0.0),
  wt(0.0),
  dz(0.0),
  uuold(0.0),
  duuold(0.0)
{}

template <class Scalar, class LocalOrdinal>
Linear2NodeFEBasis<Scalar, LocalOrdinal>::~Linear2NodeFEBasis(){}

// Calculates the values of z and u at the specified Gauss point
template<class Scalar, class LocalOrdinal>
void Linear2NodeFEBasis<Scalar, LocalOrdinal>::
compute_basis(LocalOrdinal gp, Scalar* z, Scalar* u, Scalar* uold) {
  if (gp==0) {eta=-1.0/sqrt(3.0); wt=1.0;}
  if (gp==1) {eta=1.0/sqrt(3.0); wt=1.0;}

  // Calculate basis function and derivatives at Gauss point
  phi[0]=(1.0-eta)/2.0;
  phi[1]=(1.0+eta)/2.0;
  dphide[0]=-0.5;
  dphide[1]=0.5;

  // Caculate function and derivative approximations at GP.
  dz=0.5*(z[1]-z[0]);
  zz=0.0;
  uu=0.0;
  duu=0.0;
  uuold=0.0;
  duuold=0.0;
  for (LocalOrdinal i=0; i < 2; i++) {
    zz += z[i] * phi[i];
    uu += u[i] * phi[i];
    duu += u[i] * dphide[i];
    if (uold) {
      uuold += uold[i] * phi[i];
      duuold += uold[i] * dphide[i];
    }
  }
}

#endif // FORTRILINOS_TPETRA_MODELEVALUATOR_1DFEM_DEF_HPP
