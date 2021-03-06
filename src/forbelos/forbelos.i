/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module forbelos

%include "fortrilinos_copyright.i"
%include "forerror/extern_forerror.i"
// TODO?: %import "forteuchos/forteuchos.i"
%include <std_string.i>

// Configuration
#define BELOS_DEPRECATED

// All enums should be prefaced with Belos
%rename("Belos%s", %$isenumitem) "";
%rename("Belos%s", %$isenum)     "";

// Use link-time wrappers for enums with C-only values
%nofortranconst Belos::MsgType;
%nofortranconst Belos::StatusType;
%nofortranconst Belos::ResetType;

%include "Belos_Types.i"
