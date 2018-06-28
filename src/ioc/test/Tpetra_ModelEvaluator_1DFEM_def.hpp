#ifndef TPETRA_MODELEVALUATOR_1DFEM_DEF_HPP
#define TPETRA_MODELEVALUATOR_1DFEM_DEF_HPP

// Kokkos support
#include "Kokkos_Core.hpp"

// Nonmember constuctors

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>>
tpetraModelEvaluator1DFEM(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                     const Tpetra::global_size_t num_global_elems,
                     const Scalar z_min,
                     const Scalar z_max)
{
  return Teuchos::rcp(new TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>(comm,num_global_elems,z_min,z_max));
}

// Constructor

template<class Scalar, class LO, class GO, class Node>
TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>::
TpetraModelEvaluator1DFEM(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                     const Tpetra::global_size_t num_global_elems,
                     const Scalar z_min,
                     const Scalar z_max) :
  ThyraToTpetraModelEvaluator<Scalar,LO,GO,Node>(),
  comm_(comm)
{

  TEUCHOS_ASSERT(nonnull(comm_));

  const Tpetra::global_size_t num_nodes = num_global_elems + 1;

  // owned space
  GO idx_base = 0;
  x_owned_map_ = Teuchos::rcp(new const tpetra_map(num_nodes, idx_base, comm_));

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
    x_ghosted_map_ = Teuchos::rcp(new tpetra_map(invalid, node_gids, idx_base, comm_));
  }

  importer_ = Teuchos::rcp(new Tpetra::Import<LO,GO,Node>(x_owned_map_, x_ghosted_map_));

  // residual space
  f_owned_map_ = x_owned_map_;

  // Initialize the graph for W CrsMatrix object
  graph_ = create_graph(x_owned_map_, x_ghosted_map_);

  // Create the nodal coordinates
  node_coords_ = create_mesh(x_owned_map_, z_min, z_max, num_global_elems);

  // Timers
  resid_timer_ = Teuchos::TimeMonitor::getNewCounter("Model Evaluator: Residual Evaluation");
  jac_timer_ = Teuchos::TimeMonitor::getNewCounter("Model Evaluator: Jacobian Evaluation");

  this->setup(x_owned_map_, Teuchos::null);
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Tpetra::CrsGraph<LO, GO, Node>>
TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>::
create_graph(const Teuchos::RCP<const tpetra_map>& owned_map,
             const Teuchos::RCP<const tpetra_map>& ghosted_map)
{

  // Create the shell for the graph
  Teuchos::RCP<tpetra_graph> graph =
    Teuchos::rcp(new tpetra_graph(owned_map, ghosted_map, 5));

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

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<Tpetra::MultiVector<Scalar,LO,GO,Node>>
TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>::
create_mesh(const Teuchos::RCP<const tpetra_map>& owned_map,
            const Scalar z_min, const Scalar z_max,
            const Tpetra::global_size_t num_elems)
{
  size_t num_local_nodes = owned_map->getNodeNumElements();
  GO min_GID = owned_map->getMinGlobalIndex();
  Scalar dz = (z_max - z_min)/static_cast<Scalar>(num_elems);

  Teuchos::RCP<tpetra_multi_vec> coords = Teuchos::rcp(new tpetra_multi_vec(owned_map, 1));
  for (LO node=0; node<static_cast<LO>(num_local_nodes); node++) {
    coords->replaceLocalValue(node, 0, z_min + dz * static_cast<Scalar>(min_GID+node));
  }
  return coords;
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<Tpetra::Operator<Scalar,LO,GO,Node>>
TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>::
create_tpetra_op() const
{
  Teuchos::RCP<tpetra_op> op = Teuchos::rcp(new tpetra_matrix(graph_));
  return op;
}

template<class Scalar, class LO, class GO, class Node>
void TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>::
update_x(const Teuchos::RCP<const tpetra_multi_vec>& x) const
{
  typedef Kokkos::HostSpace host_space;

  // Create ghosted objects
  if (is_null(u_ptr_))
    u_ptr_ = Teuchos::rcp(new tpetra_multi_vec(x_ghosted_map_, 1));
  u_ptr_->template modify<host_space>();
  u_ptr_->doImport(*x, *importer_, Tpetra::REPLACE);

  if (is_null(x_ptr_)) {
    x_ptr_ = Teuchos::rcp(new tpetra_multi_vec(x_ghosted_map_, 1));
    x_ptr_->template modify<host_space>();
    x_ptr_->doImport(*node_coords_, *importer_, Tpetra::INSERT);
  }

  x_ptr_->template sync<host_space>();
  u_ptr_->template sync<host_space>();
}

template<class Scalar, class LO, class GO, class Node>
void TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>::
evaluate_residual(Teuchos::RCP<tpetra_multi_vec>& f) const
{

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
    Scalar xx[2];
    xx[0] = x(ne, 0);
    xx[1] = x(ne+1, 0);

    Scalar uu[2];
    uu[0] = u(ne, 0);
    uu[1] = u(ne+1, 0);

    Linear2NodeFEBasis<Scalar, LO> basis;

    // Loop Over Gauss Points
    for(LO gp=0; gp < 2; gp++) {

      // Calculate the basis function at the gauss point
      basis.compute_basis(gp, xx, uu);

      // Loop over Nodes in Element
      for (LO i=0; i< 2; i++) {
        LO local_row =
          x_owned_map_->getLocalElement(x_ghosted_map_->getGlobalElement(ne+i));
        if (local_row != invalid) {
          Scalar value = basis.wt * basis.dz *
            (basis.uu * basis.uu * basis.phi[i] + (basis.duu * basis.dphide[i])/(basis.dz * basis.dz));
          f->sumIntoLocalValue(local_row, 0, value);
        }
      }
    }

    // Correct for Dirichlet BCs
    if ((my_rank == 0) && (ne == 0)) {
      Scalar value = u(0, 0) - 1.0;
      f->replaceLocalValue(0, 0, value);
    }

  }

}

template<class Scalar, class LO, class GO, class Node>
void TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>::
evaluate_jacobian(Teuchos::RCP<tpetra_op>& op) const
{
  Teuchos::RCP<tpetra_matrix> J = Teuchos::rcp_dynamic_cast<tpetra_matrix>(op);
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
    Scalar xx[2];
    xx[0] = x(ne, 0);
    xx[1] = x(ne+1, 0);

    Scalar uu[2];
    uu[0] = u(ne, 0);
    uu[1] = u(ne+1, 0);

    Linear2NodeFEBasis<Scalar, LO> basis;

    // Loop Over Gauss Points
    for(LO gp=0; gp < 2; gp++) {

      // Calculate the basis function at the gauss point
      basis.compute_basis(gp, xx, uu);

      // Loop over Nodes in Element
      for (LO i=0; i< 2; i++) {
        LO local_row =
          x_owned_map_->getLocalElement(x_ghosted_map_->getGlobalElement(ne+i));
        if (local_row != invalid) {
          // Loop over trial functions
          for (LO j=0; j<2; ++j) {
            const LO local_col = ne + j;
            Scalar value = basis.wt * basis.dz
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
      Scalar value = 1.0;
      J->replaceLocalValues(row, 1, &value, &column);
      column = 1;
      value = 0.0;
      J->replaceLocalValues(row, 1, &value, &column);
    }
  }

  J->fillComplete();
}

template<class Scalar, class LO, class GO, class Node>
void TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>::
evaluate_preconditioner(Teuchos::RCP<tpetra_op>& op) const
{

  Teuchos::RCP<tpetra_matrix> M_inv = Teuchos::rcp_dynamic_cast<tpetra_matrix>(op);
  TEUCHOS_ASSERT(nonnull(M_inv));
  M_inv->resumeFill();

  const LO invalid = Tpetra::Details::OrdinalTraits<LO>::invalid();
  int my_rank = comm_->getRank();
  auto num_my_elems = x_ghosted_map_->getNodeNumElements()-1;
  auto row_map = M_inv->getRowMap();
  auto col_map = M_inv->getColMap();

  if (is_null(J_diagonal_))
    J_diagonal_ = Teuchos::rcp(new tpetra_vec(x_owned_map_));
  M_inv->setAllToScalar(0.0);
  J_diagonal_->putScalar(0.0);

  Teuchos::TimeMonitor timer(*jac_timer_);

  // Loop Over # of Finite Elements on Processor
  auto x = x_ptr_->template getLocalView<Kokkos::HostSpace>();
  auto u = u_ptr_->template getLocalView<Kokkos::HostSpace>();

  for (LO ne=0; ne<static_cast<LO>(num_my_elems); ne++) {

    // Get the solution and coordinates at the nodes
    Scalar xx[2];
    xx[0] = x(ne, 0);
    xx[1] = x(ne+1, 0);

    Scalar uu[2];
    uu[0] = u(ne, 0);
    uu[1] = u(ne+1, 0);

    Linear2NodeFEBasis<Scalar, LO> basis;

    // Loop Over Gauss Points
    for(LO gp=0; gp < 2; gp++) {

      // Calculate the basis function at the gauss point
      basis.compute_basis(gp, xx, uu);

      // Loop over Nodes in Element
      for (LO i=0; i< 2; i++) {
        LO local_row =
          x_owned_map_->getLocalElement(x_ghosted_map_->getGlobalElement(ne+i));
        if (local_row != invalid) {
          // Loop over trial functions
          for (LO j=0; j<2; ++j) {
            const LO local_col = ne + j;
            if (row_map->getGlobalElement(local_row) == col_map->getGlobalElement(local_col)) {
              Scalar value = basis.wt * basis.dz
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
      Scalar value = 1.0;
      M_inv->replaceLocalValues(row, 1, &value, &column);
    }
  }

  // Invert the Jacobian diagonal for the preconditioner
  // For some reason the matrix must be fill complete before calling rightScale
  M_inv->fillComplete();
  tpetra_vec& diag = *J_diagonal_;
  M_inv->getLocalDiagCopy(diag);
  diag.reciprocal(diag);
  M_inv->rightScale(diag);
  M_inv->rightScale(diag);
}

// Finite Element Basis Object
template <class Scalar, class LO>
Linear2NodeFEBasis<Scalar, LO>::Linear2NodeFEBasis():
  uu(0.0),
  zz(0.0),
  duu(0.0),
  eta(0.0),
  wt(0.0),
  dz(0.0),
  uuold(0.0),
  duuold(0.0)
{}

template <class Scalar, class LO>
Linear2NodeFEBasis<Scalar, LO>::~Linear2NodeFEBasis(){}

// Calculates the values of z and u at the specified Gauss point
template<class Scalar, class LO>
void Linear2NodeFEBasis<Scalar, LO>::
compute_basis(LO gp, Scalar* z, Scalar* u, Scalar* uold) {
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
  for (LO i=0; i < 2; i++) {
    zz += z[i] * phi[i];
    uu += u[i] * phi[i];
    duu += u[i] * dphide[i];
    if (uold) {
      uuold += uold[i] * phi[i];
      duuold += uold[i] * dphide[i];
    }
  }
}

#endif // TPETRA_MODELEVALUATOR_1DFEM_DEF_HPP
