#ifndef FORTRILINOS_FORTRAN_OPERATOR_HPP
#define FORTRILINOS_FORTRAN_OPERATOR_HPP

#include <Teuchos_RCP.hpp>

#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_Operator.hpp>

namespace ForTrilinos {

  typedef void (*OperatorCallback)(std::pair<const double*, size_t> x, std::pair<double*, size_t> y);

#ifndef SWIG
  class FortranOperator : public Tpetra::Operator<double, int, long long, Kokkos::Compat::KokkosSerialWrapperNode> {
  private:
    typedef double                                  SC;
    typedef int                                     LO;
    typedef long long                               GO;
    typedef Kokkos::Compat::KokkosSerialWrapperNode NO;

  public:
    //! The type of the entries of the input and output multivectors.
    typedef SC scalar_type;

    //! The local index type.
    typedef LO local_ordinal_type;

    //! The global index type.
    typedef GO global_ordinal_type;

    //! The Kokkos Node type.
    typedef NO node_type;

  private:
    typedef Tpetra::Map<LO,GO,NO>               Map;
    typedef Tpetra::MultiVector<SC,LO,GO,NO>    MultiVector;

  public:
    FortranOperator() = delete;

    //! Inversion-of-Control constructor
    FortranOperator(OperatorCallback funcptr,
                    const Teuchos::RCP<const Map>& domainMap,
                    const Teuchos::RCP<const Map>& rangeMap)
        : funcptr_(funcptr)
        , domainMap_(domainMap)
        , rangeMap_(rangeMap)
    { }

    ~FortranOperator() = default;

    //! The Map associated with the domain of this operator, which must be compatible with X.getMap().
    Teuchos::RCP<const Map> getDomainMap() const {
      return domainMap_;
    }

    //! The Map associated with the range of this operator, which must be compatible with Y.getMap().
    Teuchos::RCP<const Map> getRangeMap() const {
      return rangeMap_;
    }

    //! \brief Computes the operator-multivector application.
    void apply(const MultiVector& X, MultiVector& Y,
               Teuchos::ETransp mode = Teuchos::NO_TRANS,
               SC alpha = Teuchos::ScalarTraits<SC>::one(),
               SC beta  = Teuchos::ScalarTraits<SC>::one()) const;

    /// \brief Whether this operator supports applying the transpose or conjugate transpose.
    bool hasTransposeApply() const {
      return false;
    }

  private:
    OperatorCallback funcptr_;

    Teuchos::RCP<const Map> domainMap_;
    Teuchos::RCP<const Map> rangeMap_;
  };
#endif

}

#endif // FORTRILINOS_FORTRAN_OPERATOR_HPP
