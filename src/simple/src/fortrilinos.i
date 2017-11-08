%module fortrilinos

// Load Teuchos definitions
%import <forteuchos.i>
%import <fortpetra.i>

%include "ForTrilinosSimpleInterface_config.hpp"

// Generate wrappers
%include "fortran_operator.i"
%include "solver_handle.i"
%include "eigen_handle.i"
