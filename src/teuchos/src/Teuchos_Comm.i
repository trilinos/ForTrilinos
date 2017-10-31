//---------------------------------*-SWIG-*----------------------------------//
/*!
 * \file   parameterlist/Teuchos_ParameterList.i
 * \author Seth R Johnson
 * \date   Tue Dec 06 17:59:15 2016
 * \note   Copyright (c) 2016 Oak Ridge National Laboratory, UT-Battelle, LLC.
 */
//---------------------------------------------------------------------------//
%{
#include "Teuchos_Comm.hpp"
#ifdef HAVE_MPI
# include "Teuchos_DefaultMpiComm.hpp"
#else
  typedef int MPI_Comm;
#endif
#include "Teuchos_DefaultSerialComm.hpp"
%}

%include <typemaps.i>

%include <Teuchos_RCP.i>

typedef int MPI_Comm;
%typemap(in, noblock=1) MPI_Comm %{
#ifdef HAVE_MPI
    $1 = ($1_ltype)(MPI_Comm_f2c(*(MPI_Fint *)($input)));
#else
    $1 = *$input;
#endif
%}

// Make the Comm an RCP
%teuchos_rcp(Teuchos::Comm<int>)

namespace Teuchos
{

// We provide a reduced functionality for Teuchos::comm. Ideally, we would have
// just left the constructors, and whitelisted few methods. Unfortunately, it's
// not possible to %extend a purely virtual class (and Teuchos::Comm is that).
template<class Ordinal>
class Comm
{
  public:

  %extend {
    Comm(MPI_Comm rawMpiComm) {
%#ifdef HAVE_MPI
      return static_cast<Teuchos::Comm<Ordinal>*>(new Teuchos::MpiComm<int>(rawMpiComm));
%#else
      throw std::runtime_error("MPI based constructor cannot be called when MPI is not enabled.");
%#endif
    }

    Comm() {
%#ifdef HAVE_MPI
      return static_cast<Teuchos::Comm<Ordinal>*>(new Teuchos::MpiComm<int>(MPI_COMM_WORLD));
%#else
      return static_cast<Teuchos::Comm<Ordinal>*>(new Teuchos::SerialComm<int>());
%#endif
    }

    int getRank() const {
      return $self->getRank();
    }
    int getSize() const {
      return $self->getSize();
    }
    void barrier() const {
      $self->barrier();
    }
  } // %extend
};

}

%template(TeuchosComm) Teuchos::Comm<int>;
