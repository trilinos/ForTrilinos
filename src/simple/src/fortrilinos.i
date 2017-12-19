/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module fortrilinos

%include "copyright.i"

%import <forerror.i>
%import <forteuchos.i>
%import <fortpetra.i>

%include "ForTrilinosSimpleInterface_config.hpp"

// Generate wrappers
%include "fortran_operator.i"
%include "solver_handle.i"
%include "eigen_handle.i"
