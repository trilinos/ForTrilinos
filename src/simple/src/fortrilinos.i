%module fortrilinos

// Load Teuchos definitions
%import "forteuchos.i"

%include "ForTrilinosSimpleInterface_config.hpp"

// MPI SUPPORT
// TODO: move to teuchos?
typedef int MPI_Comm;

%typemap(in, noblock=1) MPI_Comm %{
    $1 = ($1_ltype)(MPI_Comm_f2c(*(MPI_Fint *)($input)));
%}

// Generate wrappers
%include "trilinos_handle.i"

