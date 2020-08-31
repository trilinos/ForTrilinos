/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#include "model_evaluator.hpp"

#include <Teuchos_AbstractFactoryStd.hpp>

// Thyra support
#include <Thyra_DefaultSpmdVectorSpace.hpp>
#include <Thyra_DefaultSerialDenseLinearOpWithSolveFactory.hpp>
#include <Thyra_DetachedMultiVectorView.hpp>
#include <Thyra_DetachedVectorView.hpp>
#include <Thyra_MultiVectorStdOps.hpp>
#include <Thyra_VectorStdOps.hpp>
#include <Thyra_PreconditionerBase.hpp>

#include <BelosTypes.hpp>
#include <Stratimikos_DefaultLinearSolverBuilder.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_Ifpack2PreconditionerFactory.hpp>

// Tpetra support
#include <Thyra_TpetraThyraWrappers.hpp>

// Kokkos support
#include <Kokkos_Core.hpp>

namespace ForTrilinos {

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ModelEvaluator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  setup(Teuchos::RCP<Teuchos::ParameterList>& plist) {
    TEUCHOS_TEST_FOR_EXCEPTION(this->get_x_map().is_null(), std::logic_error,
        "get_x_map() must return a nonnull map when setup is called!");

    {
      auto p = Teuchos::rcpFromRef(plist->sublist("InOut Arg Settings"));
      this->setup_in_out_args(p, this->get_x_map(), this->get_f_map());
    }

    {
      auto p = Teuchos::rcpFromRef(plist->sublist("Linear Solver Settings"));
      this->setup_linear_solver(p);
    }
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ModelEvaluator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  setup_in_out_args(Teuchos::RCP<Teuchos::ParameterList>& plist,
                    const Teuchos::RCP<const Map>& x_map,
                    const Teuchos::RCP<const Map>& f_map) {
    TEUCHOS_ASSERT(nonnull(x_map));

    typedef ::Thyra::ModelEvaluatorBase MEB;
    typedef Teuchos::ScalarTraits<SC> ST;

    // Vector spaces
    TEUCHOS_ASSERT(!x_map.is_null());
    x_space_ = ::Thyra::createVectorSpace<SC,LO,GO,NO>(x_map);
    if (f_map.is_null())
      f_space_ = x_space_;
    else
      f_space_ = ::Thyra::createVectorSpace<SC,LO,GO,NO>(f_map);

    x0_ = ::Thyra::createMember(x_space_);
    V_S(x0_.ptr(), ST::zero());

    MEB::InArgsSetup<SC> in_args;
    in_args.setModelEvalDescription(this->description());
    in_args.setSupports(MEB::IN_ARG_x);
    if (plist->isParameter("Np")) {
      in_args.set_Np(plist->get<int>("Np"));
    }

    // TODO: Can call in_args.set_Np based on value in parameter list?
    prototype_in_args_ = in_args;

    MEB::OutArgsSetup<SC> out_args;
    out_args.setModelEvalDescription(this->description());
    out_args.setSupports(MEB::OUT_ARG_f);
    out_args.setSupports(MEB::OUT_ARG_W_op);
    out_args.setSupports(MEB::OUT_ARG_W_prec);
    if (plist->isParameter("Np")) {
      out_args.set_Np_Ng(plist->get<int>("Np"), 0);
    }
    prototype_out_args_ = out_args;

    nominal_values_ = in_args;
    nominal_values_.set_x(x0_);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ModelEvaluator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  setup_linear_solver(Teuchos::RCP<Teuchos::ParameterList>& plist)
  {
    Stratimikos::DefaultLinearSolverBuilder builder;

    std::string prec = plist->get("Preconditioner Type", "None");
    if (prec == "None") {
      // Do nothing
    } else if (prec == "Ifpack2") {
      using Base = Thyra::PreconditionerFactoryBase<Scalar>;
      using Impl = Thyra::Ifpack2PreconditionerFactory<Tpetra::CrsMatrix<Scalar,LO,GO,Node>>;
      builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), "Ifpack2");
    } else {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "Preconditioner Type must be one of 'None', 'Ifpack2'")
    }
    builder.setParameterList(plist);

    lows_factory = builder.createLinearSolveStrategy("");

    set_W_factory(lows_factory);

    // Create the initial guess
    initial_guess = getNominalValues().get_x()->clone_v();
    Thyra::V_S(initial_guess.ptr(),Teuchos::ScalarTraits<Scalar>::one());

  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ModelEvaluator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  set_x0(const Teuchos::ArrayView<const SC> &x0_in) {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_ASSERT_EQUALITY(x_space_->dim(), x0_in.size());
#endif
    Thyra::DetachedVectorView<SC> x0(x0_);
    x0.sv().values()().assign(x0_in);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP<Thyra::LinearOpBase<Scalar>>
  ModelEvaluator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::create_W_op() const {
    Teuchos::RCP<Operator> W_tpetra = create_operator();
    return Thyra::tpetraLinearOp<SC, LO, GO, NO>(f_space_, x_space_, W_tpetra);
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  Teuchos::RCP< ::Thyra::PreconditionerBase<Scalar>>
  ModelEvaluator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::create_W_prec() const {
    Teuchos::RCP<Operator> W_tpetra = create_operator();

    Teuchos::RCP<ThyraOp> W_op =
        Thyra::tpetraLinearOp<SC, LO, GO, NO>(f_space_, x_space_, W_tpetra);

    Teuchos::RCP<Thyra::DefaultPreconditioner<SC>> prec =
        Teuchos::rcp(new Thyra::DefaultPreconditioner<SC>);

    prec->initializeRight(W_op);
    return prec;
  }

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  void ModelEvaluator<Scalar,LocalOrdinal,GlobalOrdinal,Node>::
  evalModelImpl(const Thyra::ModelEvaluatorBase::InArgs<SC> &in_args,
                const Thyra::ModelEvaluatorBase::OutArgs<SC> &out_args) const {
    TEUCHOS_ASSERT(nonnull(in_args.get_x()));

    const Teuchos::RCP<ThyraVector> f_out = out_args.get_f();
    const Teuchos::RCP<ThyraOp> W_out = out_args.get_W_op();
    const Teuchos::RCP<ThyraPrec> W_prec_out = out_args.get_W_prec();

    const bool fill_f = nonnull(f_out);
    const bool fill_W = nonnull(W_out);
    const bool fill_W_prec = nonnull(W_prec_out);

    typedef ::Thyra::TpetraOperatorVectorExtraction<SC,LO,GO,NO> tpetra_extract;

    Teuchos::RCP<const MultiVector> x =
        tpetra_extract::getConstTpetraMultiVector(in_args.get_x());

    if (fill_f) {
      // Get the underlying tpetra objects
      Teuchos::RCP<MultiVector> f = tpetra_extract::getTpetraMultiVector(f_out);
      evaluate_residual(x, f);
    }

    if (fill_W) {
      Teuchos::RCP<Operator> W_tpetra = tpetra_extract::getTpetraOperator(W_out);
      evaluate_jacobian(x, W_tpetra);
    }

    if (fill_W_prec) {
      Teuchos::RCP<Operator> M_tpetra = tpetra_extract::getTpetraOperator(W_prec_out->getNonconstRightPrecOp());
      evaluate_preconditioner(x, M_tpetra);
    }
  }

}
