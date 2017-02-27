%module fortrilinos

// Needed to support support std::string
%{
#include <string>
#include <algorithm>
#include <stdexcept>
%}

// Load Teuchos definitions
%import "forteuchos.i"

// MPI SUPPORT
// TODO: move to teuchos?
typedef int MPI_Comm;

%typemap(in, noblock=1) MPI_Comm %{
    $1 = ($1_ltype)(MPI_Comm_f2c(*(MPI_Fint *)($input)));
%}

%apply SWIGFUNPTR { void (*)(int n, const double* x, double* y) } ;

// Generate wrappers
%include "trilinos_handle.i"

