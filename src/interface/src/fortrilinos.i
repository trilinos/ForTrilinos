/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module fortrilinos

%include <copyright.i>
%include <extern_forerror.i>

%import <forteuchos.i>
%import <fortpetra.i>

%include "ForTrilinosInterface_config.hpp"

%inline %{
#include "Kokkos_DefaultNode.hpp"

typedef double                                  SC;
typedef int                                     LO;
typedef long long                               GO;
typedef Kokkos::Compat::KokkosSerialWrapperNode NO;
typedef char                                    Packet;
%}

// Generate wrappers
%{
#include "solver_handle.hpp"
#include "eigen_handle.hpp"
%}

%include "solver_handle.hpp"
%include "eigen_handle.hpp"

%include "nox.i"
