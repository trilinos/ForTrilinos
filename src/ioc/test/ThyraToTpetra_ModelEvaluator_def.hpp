#ifndef THYRATOTPETRA_MODELEVALUATOR_DEF_HPP
#define THYRATOTPETRA_MODELEVALUATOR_DEF_HPP

// Thyra support
#include "Thyra_DefaultSpmdVectorSpace.hpp"
#include "Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp"
#include "Thyra_DetachedMultiVectorView.hpp"
#include "Thyra_DetachedVectorView.hpp"
#include "Thyra_MultiVectorStdOps.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_PreconditionerBase.hpp"

// Tpetra support
#include "Thyra_TpetraThyraWrappers.hpp"

// Kokkos support
#include "Kokkos_Core.hpp"

// Constructor
template<class Scalar, class LO, class GO, class Node>
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::
ThyraToTpetraModelEvaluator() : show_get_invalid_arg_(false) { }

// Initializers/Accessors

template<class Scalar, class LO, class GO, class Node>
void
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::
setup(const Teuchos::RCP<const tpetra_map>& x_map,
      const Teuchos::RCP<const tpetra_map>& f_map)
{

  TEUCHOS_ASSERT(nonnull(x_map));

  typedef ::Thyra::ModelEvaluatorBase MEB;
  typedef Teuchos::ScalarTraits<scalar_type> ST;

  // Vector spaces
  x_space_ = ::Thyra::createVectorSpace<Scalar, LO, GO, Node>(x_map);
  if (f_map.is_null())
    f_space_ = x_space_;
  else
    f_space_ = ::Thyra::createVectorSpace<Scalar, LO, GO, Node>(f_map);

  x0_ = ::Thyra::createMember(x_space_);
  V_S(x0_.ptr(), ST::zero());

  MEB::InArgsSetup<Scalar> in_args;
  in_args.setModelEvalDescription(this->description());
  in_args.setSupports(MEB::IN_ARG_x);
  prototype_in_args_ = in_args;

  MEB::OutArgsSetup<Scalar> out_args;
  out_args.setModelEvalDescription(this->description());
  out_args.setSupports(MEB::OUT_ARG_f);
  out_args.setSupports(MEB::OUT_ARG_W_op);
  out_args.setSupports(MEB::OUT_ARG_W_prec);
  prototype_out_args_ = out_args;

  nominal_values_ = in_args;
  nominal_values_.set_x(x0_);
}

template<class Scalar, class LO, class GO, class Node>
void ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::
evaluate_jacobian(Teuchos::RCP<tpetra_op>& op) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      "evaluate_jacobian must be implemented by derived classes");
}

template<class Scalar, class LO, class GO, class Node>
void ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::
evaluate_preconditioner(Teuchos::RCP<tpetra_op>& op) const
{
  TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
      "evaluate_preconditioner must be implemented by derived classes");
}

template<class Scalar, class LO, class GO, class Node>
void ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::
set_x0(const Teuchos::ArrayView<const Scalar> &x0_in)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
  Thyra::DetachedVectorView<Scalar> x0(x0_);
  x0.sv().values()().assign(x0_in);
}

template<class Scalar, class LO, class GO, class Node>
void ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::
setShowGetInvalidArgs(bool show_get_invalid_arg)
{
  show_get_invalid_arg_ = show_get_invalid_arg;
}

template<class Scalar, class LO, class GO, class Node>
void ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::
set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<Scalar>>& W_factory)
{
  W_factory_ = W_factory;
}

// Public functions overridden from ModelEvaulator


template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::get_x_space() const
{
  return x_space_;
}


template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Thyra::VectorSpaceBase<Scalar>>
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::get_f_space() const
{
  return f_space_;
}


template<class Scalar, class LO, class GO, class Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::getNominalValues() const
{
  return nominal_values_;
}


template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::create_W_op() const
{
  Teuchos::RCP<tpetra_op> W_tpetra = create_tpetra_op();
  return Thyra::tpetraLinearOp<Scalar, LO, GO, Node>(f_space_, x_space_, W_tpetra);
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP< ::Thyra::PreconditionerBase<Scalar>>
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::create_W_prec() const
{
  Teuchos::RCP<tpetra_op> W_tpetra = create_tpetra_op();

  Teuchos::RCP<thyra_op> W_op =
    Thyra::tpetraLinearOp<Scalar, LO, GO, Node>(f_space_, x_space_, W_tpetra);

  Teuchos::RCP<Thyra::DefaultPreconditioner<Scalar>> prec =
    Teuchos::rcp(new Thyra::DefaultPreconditioner<Scalar>);

  prec->initializeRight(W_op);
  return prec;
}

template<class Scalar, class LO, class GO, class Node>
Teuchos::RCP<const Thyra::LinearOpWithSolveFactoryBase<Scalar>>
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::get_W_factory() const
{
  return W_factory_;
}


template<class Scalar, class LO, class GO, class Node>
Thyra::ModelEvaluatorBase::InArgs<Scalar>
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::createInArgs() const
{
  return prototype_in_args_;
}

// Private functions overridden from ModelEvaulatorDefaultBase


template<class Scalar, class LO, class GO, class Node>
Thyra::ModelEvaluatorBase::OutArgs<Scalar>
ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::createOutArgsImpl() const
{
  return prototype_out_args_;
}

template<class Scalar, class LO, class GO, class Node>
void ThyraToTpetraModelEvaluator<Scalar, LO, GO, Node>::
evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<Scalar> &in_args,
              const Thyra::ModelEvaluatorBase::OutArgs<Scalar> &out_args) const
{
  TEUCHOS_ASSERT(nonnull(in_args.get_x()));

  const Teuchos::RCP<thyra_vec> f_out = out_args.get_f();
  const Teuchos::RCP<thyra_op> W_out = out_args.get_W_op();
  const Teuchos::RCP<thyra_prec> W_prec_out = out_args.get_W_prec();

  const bool fill_f = nonnull(f_out);
  const bool fill_W = nonnull(W_out);
  const bool fill_W_prec = nonnull(W_prec_out);

  typedef Tpetra::Operator<Scalar,LO,GO,Node> tpetra_op;
  typedef ::Thyra::TpetraOperatorVectorExtraction<Scalar,LO,GO,Node> tpetra_extract;

  Teuchos::RCP<const tpetra_multi_vec> x =
    tpetra_extract::getConstTpetraMultiVector(in_args.get_x());

  update_x(x);

  if (fill_f) {
    // Get the underlying tpetra objects
    Teuchos::RCP<tpetra_multi_vec> f = tpetra_extract::getTpetraMultiVector(f_out);
    evaluate_residual(f);
  }

  if (fill_W) {
    Teuchos::RCP<tpetra_op> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
    evaluate_jacobian(W_tpetra);
  }

  if (fill_W_prec) {
    Teuchos::RCP<tpetra_op> M_tpetra = tpetra_extract::getTpetraOperator(W_prec_out->getNonconstRightPrecOp());
    evaluate_preconditioner(M_tpetra);
  }

}

#endif // THYRATOTPETRA_MODELEVALUATOR_DEF_HPP
