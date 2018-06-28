#ifndef THYRATOTPETRA_MODELEVALUATOR_HPP
#define THYRATOTPETRA_MODELEVALUATOR_HPP

#include "Thyra_StateFuncModelEvaluatorBase.hpp"

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Import.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Teuchos_TimeMonitor.hpp"

/** \brief Base class for ThyraToTpetra model evaluators

 */
template<class Scalar, class LO, class GO, class Node>
class ThyraToTpetraModelEvaluator
  : public ::Thyra::StateFuncModelEvaluatorBase<Scalar>
{
public:

  // Public typedefs
  typedef Scalar scalar_type;
  typedef LO local_ordinal_type;
  typedef GO global_ordinal_type;
  typedef Node node_type;
  typedef Tpetra::Map<LO, GO, Node> tpetra_map;
  typedef Tpetra::Operator<Scalar, LO, GO, Node> tpetra_op;
  typedef Tpetra::MultiVector<Scalar, LO, GO, Node> tpetra_multi_vec;
  typedef ::Thyra::VectorBase<Scalar> thyra_vec;
  typedef ::Thyra::VectorSpaceBase<Scalar> thyra_vec_space;
  typedef ::Thyra::LinearOpBase<Scalar> thyra_op;
  typedef ::Thyra::PreconditionerBase<Scalar> thyra_prec;

  // Constructor
  ThyraToTpetraModelEvaluator();

  /** \name Initializers/Accessors */
  //@{

  /** \brief . */
  void set_x0(const Teuchos::ArrayView<const Scalar> &x0);

  /** \brief . */
  void setShowGetInvalidArgs(bool show_get_invalid_arg);

  void set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar>>& W_factory);

  //@}

  /** \name Public functions overridden from ModelEvaulator. */
  //@{

  /** \brief . */
  Teuchos::RCP<const thyra_vec_space> get_x_space() const;

  /** \brief . */
  Teuchos::RCP<const thyra_vec_space> get_f_space() const;

  /** \brief . */
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> getNominalValues() const;

  /** \brief . */
  Teuchos::RCP<thyra_op> create_W_op() const;

  /** \brief . */
  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar>> get_W_factory() const;

  /** \brief . */
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> createInArgs() const;

  /** \brief . */
  Teuchos::RCP<thyra_prec> create_W_prec() const;
  //@}

  void
  setup(const Teuchos::RCP<const tpetra_map>& x_map,
        const Teuchos::RCP<const tpetra_map>& f_map);

  virtual void
  evaluate_residual(Teuchos::RCP<tpetra_multi_vec>& f) const = 0;

  virtual void
  evaluate_jacobian(Teuchos::RCP<tpetra_op>& J) const;

  virtual void
  evaluate_preconditioner(Teuchos::RCP<tpetra_op>& J) const;

  virtual void
  update_x(const Teuchos::RCP<const tpetra_multi_vec>& x) const = 0;

  virtual Teuchos::RCP<tpetra_op>
  create_tpetra_op() const = 0;

private:

  /** \name Private functions overridden from ModelEvaulatorDefaultBase. */
  //@{

  /** \brief . */
  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> createOutArgsImpl() const;

  /** \brief . */
  void evalModelImpl(
    const ::Thyra::ModelEvaluatorBase::InArgs<Scalar> &in_args,
    const ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> &out_args
    ) const;

  //@}

private: // data members

  Teuchos::RCP<const thyra_vec_space> x_space_;
  Teuchos::RCP<const tpetra_map>   x_map_;

  Teuchos::RCP<const thyra_vec_space> f_space_;
  Teuchos::RCP<const tpetra_map>   f_map_;

  Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar>> W_factory_;

  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> nominal_values_;
  Teuchos::RCP<thyra_vec> x0_;
  bool show_get_invalid_arg_;
  ::Thyra::ModelEvaluatorBase::InArgs<Scalar> prototype_in_args_;
  ::Thyra::ModelEvaluatorBase::OutArgs<Scalar> prototype_out_args_;

};

//==================================================================
#include "ThyraToTpetra_ModelEvaluator_def.hpp"
//==================================================================

#endif // THYRATOTPETRA_MODELEVALUATOR_HPP
