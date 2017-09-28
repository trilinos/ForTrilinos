%module fortrilinos

// Load Teuchos definitions
%import "forteuchos.i"

%include "ForTrilinosSimpleInterface_config.hpp"

// Generate wrappers
%include "solver_handle.i"
%include "eigen_handle.i"
