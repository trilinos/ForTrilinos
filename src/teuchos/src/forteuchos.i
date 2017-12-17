/*
 * Copyright 2017, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module forteuchos

%include "copyright.i"

%include "ForTrilinosTeuchos_config.hpp"

%include <std_vector.i>
%template(VectorInt)    std::vector<int>;
%template(VectorDouble) std::vector<double>;
%template(VectorLongLong) std::vector<long long>;

// Typedefs
typedef int Teuchos_Ordinal;

%include <typemaps.i>
%fortran_view(int)
%fortran_view(long long)
%fortran_view(double)
%fortran_view(size_t)

// enum workaround
#define RENAME_ENUM(X) %rename(Teuchos##X) X;
RENAME_ENUM(ESide)
RENAME_ENUM(LEFT_SIDE)
RENAME_ENUM(RIGHT_SIDE)
RENAME_ENUM(ETransp)
RENAME_ENUM(NO_TRANS)
RENAME_ENUM(TRANS)
RENAME_ENUM(CONJ_TRANS)
RENAME_ENUM(EUplo)
RENAME_ENUM(UPPER_TRI)
RENAME_ENUM(LOWER_TRI)
RENAME_ENUM(UNDEF_TRI)
RENAME_ENUM(EDiag)
RENAME_ENUM(UNIT_DIAG)
RENAME_ENUM(NON_UNIT_DIAG)
RENAME_ENUM(EType)
RENAME_ENUM(FULL)
RENAME_ENUM(LOWER)
RENAME_ENUM(UPPER)
RENAME_ENUM(HESSENBERG)
RENAME_ENUM(SYM_BAND_L)
RENAME_ENUM(SYM_BAND_U)
RENAME_ENUM(BAND)
RENAME_ENUM(DataAccess)
RENAME_ENUM(Copy)
RENAME_ENUM(View)
#undef RENAME_ENUM

%include "Teuchos_Types.i"
%include "Teuchos_Exceptions.i"
%include "Teuchos_RCP.i"

%include "Teuchos_Comm.i"
%include "Teuchos_ParameterList.i"
%include "Teuchos_XML.i"
