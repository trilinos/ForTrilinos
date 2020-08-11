/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#ifndef FORTRILINOS_TPETRA_MODELEVALUATOR_1DFEM_HPP
#define FORTRILINOS_TPETRA_MODELEVALUATOR_1DFEM_HPP

#include <Tpetra_Map.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_CrsGraph.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "fortrilinos_hl/model_evaluator.hpp"

template<class SC, class LO, class GO, class NO>
class TpetraModelEvaluator1DFEM;

/** \brief 1D Finite Element model for nonlinear heat conduction
 *
 * The equation modeled is:

 \verbatim

   d2T
   --- - T**2 = 0
   dz2

   subject to:
      T  = 1.0 @ z = z_min
      T' = 0.0 @ z = z_max

 \endverbatim

 */
template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
class TpetraModelEvaluator1DFEM
  : public ForTrilinos::ModelEvaluator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
public:

  typedef Scalar SC;
  typedef LocalOrdinal LO;
  typedef GlobalOrdinal GO;
  typedef Node NO;

  // Public typedefs
  typedef SC scalar_type;
  typedef LO local_ordinal_type;
  typedef GO global_ordinal_type;
  typedef NO node_type;

  typedef Tpetra::CrsMatrix<SC,LO,GO,NO>    Matrix;
  typedef Tpetra::CrsGraph<LO,GO,NO>        Graph;
  typedef Tpetra::Map<LO,GO,NO>             Map;
  typedef Tpetra::MultiVector<SC,LO,GO,NO>  MultiVector;
  typedef Tpetra::Vector<SC,LO,GO,NO>       Vector;
  typedef Tpetra::Operator<SC,LO,GO,NO>     Operator;


  // Constructor
  TpetraModelEvaluator1DFEM(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                            const Tpetra::global_size_t num_global_elems,
                            const SC z_min, const SC z_max);

  /** Evaluate the FE residual */
  void evaluate_residual(const Teuchos::RCP<const MultiVector>&,
                         Teuchos::RCP<MultiVector>&) const override;

  /** Evaluate the stiffness (Jacobian) */
  void evaluate_jacobian(const Teuchos::RCP<const MultiVector>&,
                         Teuchos::RCP<Operator>&) const override;

  /** Evaluate the preconditioner */
  void evaluate_preconditioner(const Teuchos::RCP<const MultiVector>&,
                               Teuchos::RCP<Operator>&) const override;

  /** Create the Tpetra operator **/
  Teuchos::RCP<Operator> create_operator() const override;

private:

  /** Update the solution vector **/
  void update_solution_vector(const Teuchos::RCP<const MultiVector>&) const;

  /** Creates the 1D mesh */
  Teuchos::RCP<MultiVector>
  create_mesh(const Teuchos::RCP<const Map>& owned_map,
              const SC z_min, const SC z_max,
              const Tpetra::global_size_t num_elems);

  /** Allocates and returns the Jacobian matrix graph */
  Teuchos::RCP<const Graph>
  create_graph(const Teuchos::RCP<const Map>& owned_map,
               const Teuchos::RCP<const Map>& ghosted_map);

private: // data members

  const Teuchos::RCP<const Teuchos::Comm<int>>  comm_;

  Teuchos::RCP<const Graph> graph_;

  Teuchos::RCP<const Map>   x_owned_map_;
  Teuchos::RCP<const Map>   x_ghosted_map_;
  Teuchos::RCP<const Tpetra::Import<LO, GO, NO>> importer_;

  Teuchos::RCP<const Map>   f_owned_map_;

  Teuchos::RCP<MultiVector> node_coords_;

  mutable Teuchos::RCP<MultiVector> u_ptr_;
  mutable Teuchos::RCP<MultiVector> x_ptr_;
  mutable Teuchos::RCP<Vector> J_diagonal_;
  mutable Teuchos::RCP<Teuchos::Time> resid_timer_;
  mutable Teuchos::RCP<Teuchos::Time> jac_timer_;

  Teuchos::RCP<const Map> get_x_map() const override {return x_owned_map_;}
  Teuchos::RCP<const Map> get_f_map() const override {return Teuchos::null;}
};


// Finite Element Basis Object
template<class Scalar, class LocalOrdinal>
class Linear2NodeFEBasis {

 public:
  // Constructor
  Linear2NodeFEBasis();

  // Destructor
  ~Linear2NodeFEBasis();

  // Calculates the values of u and x at the specified gauss point
  void compute_basis(LocalOrdinal gp, Scalar *x, Scalar *u, Scalar *uold = 0);

 public:
  // Variables that are calculated at the gauss point
  Scalar phi[2], dphide[2];
  Scalar uu, zz, duu, eta, wt;
  Scalar dz;
  // These are only needed for transient
  Scalar uuold, duuold;
};

//==================================================================
#include "Tpetra_ModelEvaluator_1DFEM_def.hpp"
//==================================================================

#endif // FORTRILINOS_TPETRA_MODELEVALUATOR_1DFEM_HPP
