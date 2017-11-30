/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module forbelos

%include "copyright.i"

%include <std_string.i>

%include "ForTrilinosBelos_config.hpp"

%ignore Belos::toString;

#define BELOS_DEPRECATED
%include "Belos_Types.i"
