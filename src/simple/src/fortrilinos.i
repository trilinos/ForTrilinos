%module fortrilinos

// Load Teuchos definitions
%import "forteuchos.i"

%include "ForTrilinosSimpleInterface_config.hpp"

typedef int MPI_Comm;

%typemap(in, noblock=1) MPI_Comm %{
#ifdef HAVE_MPI
    $1 = ($1_ltype)(MPI_Comm_f2c(*(MPI_Fint *)($input)));
#else
    $1 = *$input;
#endif
%}

// Generate wrappers
%include "trilinos_handle.i"

