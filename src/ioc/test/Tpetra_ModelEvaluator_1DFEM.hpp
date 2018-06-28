#ifndef TPETRA_MODELEVALUATOR_1DFEM_HPP
#define TPETRA_MODELEVALUATOR_1DFEM_HPP

#include "Tpetra_Map.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_TimeMonitor.hpp"

#include "ThyraToTpetra_ModelEvaluator.hpp"

template<class Scalar, class LO, class GO, class Node>
class TpetraModelEvaluator1DFEM;

/** \brief Nonmember constuctor.
 *
 * \relates TpetraModelEvaluator1DFEM
 */
template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<TpetraModelEvaluator1DFEM<Scalar, LO, GO, Node>>
tpetraModelEvaluator1DFEM(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                          const Tpetra::global_size_t num_global_elems,
                          const Scalar z_min, const Scalar z_max);


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
template<class Scalar, class LO, class GO, class Node>
class TpetraModelEvaluator1DFEM
  : public ThyraToTpetraModelEvaluator<Scalar,LO,GO,Node>
{
public:

  // Public typedefs
  typedef Scalar scalar_type;
  typedef LO local_ordinal_type;
  typedef GO global_ordinal_type;
  typedef Node node_type;
  typedef Tpetra::Map<LO, GO, Node> tpetra_map;
  typedef Tpetra::CrsGraph<LO, GO, Node> tpetra_graph;
  typedef Tpetra::Operator<Scalar, LO, GO, Node> tpetra_op;
  typedef Tpetra::CrsMatrix<Scalar, LO, GO, Node> tpetra_matrix;
  typedef Tpetra::Vector<Scalar, LO, GO, Node> tpetra_vec;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> tpetra_multi_vec;

  // Constructor
  TpetraModelEvaluator1DFEM(const Teuchos::RCP<const Teuchos::Comm<int>>& comm,
                       const Tpetra::global_size_t num_global_elems,
                       const Scalar z_min, const Scalar z_max);

  /** Evaluate the FE residual */
  void evaluate_residual(Teuchos::RCP<tpetra_multi_vec>&) const;

  /** Evaluate the stiffness (Jacobian) */
  void evaluate_jacobian(Teuchos::RCP<tpetra_op>&) const;

  /** Evaluate the preconditioner */
  void evaluate_preconditioner(Teuchos::RCP<tpetra_op>&) const;

  void update_x(const Teuchos::RCP<const tpetra_multi_vec>& x) const;

  /** Create the Tpetra operator **/
  Teuchos::RCP<tpetra_op> create_tpetra_op() const;

private:

  /** Creates the 1D mesh */
  virtual Teuchos::RCP<tpetra_multi_vec>
  create_mesh(const Teuchos::RCP<const tpetra_map>& owned_map,
              const Scalar z_min, const Scalar z_max,
              const Tpetra::global_size_t num_elems);

  /** Allocates and returns the Jacobian matrix graph */
  virtual Teuchos::RCP<const tpetra_graph>
  create_graph(const Teuchos::RCP<const tpetra_map>& owned_map,
               const Teuchos::RCP<const tpetra_map>& ghosted_map);

private: // data members

  const Teuchos::RCP<const Teuchos::Comm<int>>  comm_;

  Teuchos::RCP<const tpetra_graph> graph_;

  Teuchos::RCP<const tpetra_map>   x_owned_map_;
  Teuchos::RCP<const tpetra_map>   x_ghosted_map_;
  Teuchos::RCP<const Tpetra::Import<LO, GO, Node>> importer_;

  Teuchos::RCP<const tpetra_map>   f_owned_map_;

  Teuchos::RCP<tpetra_multi_vec> node_coords_;

  mutable Teuchos::RCP<tpetra_multi_vec> u_ptr_;
  mutable Teuchos::RCP<tpetra_multi_vec> x_ptr_;
  mutable Teuchos::RCP<tpetra_vec> J_diagonal_;
  mutable Teuchos::RCP<Teuchos::Time> resid_timer_;
  mutable Teuchos::RCP<Teuchos::Time> jac_timer_;
};


// Finite Element Basis Object
template<class Scalar, class LO>
class Linear2NodeFEBasis {

 public:
  // Constructor
  Linear2NodeFEBasis();

  // Destructor
  ~Linear2NodeFEBasis();

  // Calculates the values of u and x at the specified gauss point
  void compute_basis(LO gp, Scalar *x, Scalar *u, Scalar *uold = 0);

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

#endif // TPETRA_MODELEVALUATOR_1DFEM_HPP
