%module fortrilinos

// Load Teuchos definitions
%import "forteuchos.i"

%include "ForTrilinosSimpleInterface_config.hpp"

// Mpi_Comm typemaps
%typemap(in, noblock=1) MPI_Comm %{
#ifdef HAVE_MPI
    $1 = ($1_ltype)(MPI_Comm_f2c(*(MPI_Fint *)($input)));
#else
    $1 = *$input;
#endif
%}
%typemap(ftype)  MPI_Comm %{integer(C_INT), intent(in)%}
%typemap(imtype) MPI_Comm %{integer(C_INT), intent(in)%}
%typemap(fin)    MPI_Comm %{$1_name%}

// Generate wrappers
%include "trilinos_handle.i"

