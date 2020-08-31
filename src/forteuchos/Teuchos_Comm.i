/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include <Teuchos_Comm.hpp>
#if FORTRILINOS_USE_MPI
# include "Teuchos_DefaultMpiComm.hpp"
#else
  typedef int MPI_Comm;
#endif
#include <Teuchos_DefaultSerialComm.hpp>
%}

%include <typemaps.i>

%include <Teuchos_RCP.i>

%apply int { MPI_Comm };
%typemap(ftype) MPI_Comm
   "integer"
%typemap(fin, noblock=1) MPI_Comm {
    $1 = int($input, C_INT)
}
%typemap(fout, noblock=1) MPI_Comm {
    $result = int($1)
}

%typemap(in, noblock=1) MPI_Comm {
%#if FORTRILINOS_USE_MPI
    $1 = MPI_Comm_f2c(%static_cast(*$input, MPI_Fint));
%#else
    $1 = *$input;
%#endif
}
%typemap(out, noblock=1) MPI_Comm {
%#if FORTRILINOS_USE_MPI
    $result = %static_cast(MPI_Comm_c2f($1), int);
%#else
    $result = $1;
%#endif
}

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

    // Wrap built-in methods
    int getRank() const;
    int getSize() const;
    void barrier() const;

    // Add constructors
  %extend {
    Comm(MPI_Comm rawMpiComm) {
%#if FORTRILINOS_USE_MPI
      return new Teuchos::MpiComm<Ordinal>(rawMpiComm);
%#else
      throw std::runtime_error("MPI based constructor cannot be called when MPI is not enabled.");
%#endif
    }

    Comm() {
%#if FORTRILINOS_USE_MPI
      return new Teuchos::MpiComm<Ordinal>(MPI_COMM_WORLD);
%#else
      return new Teuchos::SerialComm<Ordinal>();
%#endif
    }

    MPI_Comm getRawMpiComm() {
%#if FORTRILINOS_USE_MPI
      Teuchos::MpiComm<Ordinal>& comm = dynamic_cast<Teuchos::MpiComm<Ordinal>&>(*$self);
      return *comm.getRawMpiComm();
%#else
      throw std::runtime_error("MPI based constructor cannot be called when MPI is not enabled.");
%#endif
    }
  } // %extend
};

}

%template(TeuchosComm) Teuchos::Comm<int>;
