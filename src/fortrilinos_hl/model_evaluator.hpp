/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
#ifndef FORTRILINOS_MODEL_EVALUATOR_HPP
#define FORTRILINOS_MODEL_EVALUATOR_HPP

#include "ForTrilinos_config.h"

#include <Thyra_StateFuncModelEvaluatorBase.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Import.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Teuchos_TimeMonitor.hpp>

namespace ForTrilinos {

  template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  class ModelEvaluator : public ::Thyra::StateFuncModelEvaluatorBase<double> {
  public:
    typedef Scalar                                  SC;
    typedef LocalOrdinal                            LO;
    typedef GlobalOrdinal                           GO;
    typedef Node                                    NO;
    typedef size_t                                  global_size_t;

    typedef Tpetra::CrsMatrix<SC,LO,GO,NO>          Matrix;
    typedef Tpetra::Map<LO,GO,NO>                   Map;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>        MultiVector;
    typedef Tpetra::Operator<SC,LO,GO,NO>           Operator;

    typedef ::Thyra::VectorBase<SC>                 ThyraVector;
    typedef ::Thyra::VectorSpaceBase<SC>            ThyraVectorSpace;
    typedef ::Thyra::LinearOpBase<SC>               ThyraOp;
    typedef ::Thyra::PreconditionerBase<SC>         ThyraPrec;

  public:
    virtual void evaluate_residual(const Teuchos::RCP<const MultiVector>& x,
                                   Teuchos::RCP<MultiVector>& f) const = 0;

    virtual void evaluate_jacobian(const Teuchos::RCP<const MultiVector>& x,
                                   Teuchos::RCP<Operator>& J) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "evaluate_jacobian must be implemented by derived classes");
    }

    virtual void evaluate_preconditioner(const Teuchos::RCP<const MultiVector>& x,
                                         Teuchos::RCP<Operator>& M) const
    {
      TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument,
        "evaluate_preconditioner must be implemented by derived classes");
    }

    virtual Teuchos::RCP<const Map> get_x_map() const = 0;
    virtual Teuchos::RCP<const Map> get_f_map() const = 0;
    virtual Teuchos::RCP<Operator> create_operator() const = 0;

#ifndef SWIG
  public:
    Teuchos::RCP<Thyra::VectorBase<Scalar>> initial_guess;
    Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar>> lows_factory;

    // FIXME: These should not be public
    void set_W_factory(const Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<SC>>& W_factory) {
      W_factory_ = W_factory;
    }
    Teuchos::RCP<ThyraOp> create_W_op() const;

    ::Thyra::ModelEvaluatorBase::InArgs<SC> getNominalValues() const {
      return nominal_values_;
    }
#endif

    void setup(Teuchos::RCP<Teuchos::ParameterList>& plist);

  protected:
    ModelEvaluator()
        : show_get_invalid_arg_(false)
    { }

    void setShowGetInvalidArgs(bool show_get_invalid_arg) {
      show_get_invalid_arg_ = show_get_invalid_arg;
    }

    Teuchos::RCP<const ThyraVectorSpace> get_x_space() const {
      return x_space_;
    }

    Teuchos::RCP<const ThyraVectorSpace> get_f_space() const {
      return f_space_;
    }

    Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<SC>> get_W_factory() const {
      return W_factory_;
    }

    ::Thyra::ModelEvaluatorBase::InArgs<SC> createInArgs() const {
      return prototype_in_args_;
    }

    Teuchos::RCP<ThyraPrec> create_W_prec() const;


  private:

    void set_x0(const Teuchos::ArrayView<const SC> &x0);

    // Private functions overridden from ModelEvaulatorDefaultBase. */
    ::Thyra::ModelEvaluatorBase::OutArgs<SC> createOutArgsImpl() const {
      return prototype_out_args_;
    }

    void evalModelImpl(const ::Thyra::ModelEvaluatorBase::InArgs<SC> &in_args,
                       const ::Thyra::ModelEvaluatorBase::OutArgs<SC> &out_args) const;

  private:

    void setup_in_out_args(Teuchos::RCP<Teuchos::ParameterList>& plist,
                           const Teuchos::RCP<const Map>& x_map,
                           const Teuchos::RCP<const Map>& f_map);

    void setup_linear_solver(Teuchos::RCP<Teuchos::ParameterList>& plist);

  private:

    Teuchos::RCP<const ThyraVectorSpace> x_space_;
    Teuchos::RCP<const Map>              x_map_;
    Teuchos::RCP<const ThyraVectorSpace> f_space_;
    Teuchos::RCP<const Map>              f_map_;

    Teuchos::RCP<const ::Thyra::LinearOpWithSolveFactoryBase<SC>> W_factory_;

    ::Thyra::ModelEvaluatorBase::InArgs<SC> nominal_values_;
    Teuchos::RCP<ThyraVector> x0_;
    bool show_get_invalid_arg_;
    ::Thyra::ModelEvaluatorBase::InArgs<SC> prototype_in_args_;
    ::Thyra::ModelEvaluatorBase::OutArgs<SC> prototype_out_args_;
  };
}

#ifndef SWIG
#include "model_evaluator_def.hpp"
#endif

#endif // FORTRILINOS_MODEL_EVALUATOR_HPP
