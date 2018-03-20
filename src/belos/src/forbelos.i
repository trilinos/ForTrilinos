/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module forbelos

%include <copyright.i>
%include <extern_forerror.i>
%import <forerror.i>
// TODO?: %import <forteuchos.i>
%include <std_string.i>

// Configuration
#define BELOS_DEPRECATED

// All enums should be prefaced with Belos
%rename("Belos%s", %$isenumitem) "";
%rename("Belos%s", %$isenum)     "";

%include "Belos_Types.i"
