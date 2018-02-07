/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
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
%typemap(in, noblock=1) MPI_Comm {
%#ifdef HAVE_MPI
    $1 = MPI_Comm_f2c(%static_cast(*$input, MPI_Fint));
%#else
    $1 = *$input;
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

    }
  } // %extend
};

}

%template(TeuchosComm) Teuchos::Comm<int>;
