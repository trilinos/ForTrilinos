/* Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include <BelosTypes.hpp>
%}

%ignore Belos::toString;
%ignore Belos::BelosError;

%include "BelosTypes.hpp"
