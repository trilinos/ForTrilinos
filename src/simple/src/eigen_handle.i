/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%{
#include "eigen_handle.hpp"
%}

%ignore setup_matrix(Teuchos::RCP<Matrix>);
%ignore setup_matrix_rhs(Teuchos::RCP<Matrix>);

%include "Teuchos_Comm.i"
%include "eigen_handle.hpp"
