/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module fortrilinos_hl

%include "fortrilinos_copyright.i"
%include "forerror/extern_forerror.i"

%import "forteuchos/forteuchos.i"
%import "fortpetra/fortpetra.i"

%include "ForTrilinos_config.h"

%{
#include <Kokkos_DefaultNode.hpp>
#include "fortpetra/ForTrilinos_DefaultNodeType.hpp"
%}
%inline %{
typedef double                                  SC;
typedef int                                     LO;
typedef long long                               GO;
typedef ForTrilinos::DefaultNodeType            NO;
typedef char                                    Packet;
%}

// Generate wrappers
%{
#include "fortrilinos_hl/solver_handle.hpp"
#include "fortrilinos_hl/eigen_handle.hpp"
%}

%include "solver_handle.hpp"
%include "eigen_handle.hpp"

%include "nox.i"
