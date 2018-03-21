/*
 * Copyright 2017-2018, UT-Battelle, LLC
 *
 * SPDX-License-Identifier: BSD-3-Clause
 * License-Filename: LICENSE
 */
%module forbelos

%include <copyright.i>
%include <extern_forerror.i>

// Convert all std::string references/values to and from Fortran strings
%ignore std::string;
%include <std_string.i>
%apply std::string NATIVE { std::string  };
%apply const std::string& NATIVE { const std::string& };

%include "ForTrilinosBelos_config.hpp"

%ignore Belos::toString;

// enum workaround
// TODO: you can probably simplify this by using an expression like:
// %rename("Belos%s", %$isenumitem, regextarget=1, fullname=1) "Belos::.*";
#define RENAME_ENUM(X) %rename(Belos##X) X;
RENAME_ENUM(ETrans)
RENAME_ENUM(NOTRANS)
RENAME_ENUM(TRANS)
RENAME_ENUM(CONJTRANS)
RENAME_ENUM(NormType)
RENAME_ENUM(OneNorm)
RENAME_ENUM(TwoNorm)
RENAME_ENUM(InfNorm)
RENAME_ENUM(ScaleType)
RENAME_ENUM(NormOfRHS)
RENAME_ENUM(NormOfInitRes)
RENAME_ENUM(NormOfPrecInitRes)
RENAME_ENUM(None)
RENAME_ENUM(UserProvided)
RENAME_ENUM(NormOfFullInitRes)
RENAME_ENUM(NormOfFullPrecInitRes)
RENAME_ENUM(NormOfFullScaledInitRes)
RENAME_ENUM(NormOfFullScaledPrecInitRes)
RENAME_ENUM(OutputType)
RENAME_ENUM(General)
RENAME_ENUM(Brief)
RENAME_ENUM(User)
RENAME_ENUM(ReturnType)
RENAME_ENUM(Converged)
RENAME_ENUM(Unconverged)
RENAME_ENUM(StatusType)
RENAME_ENUM(Passed)
RENAME_ENUM(Failed)
RENAME_ENUM(Undefined)
RENAME_ENUM(ResetType)
RENAME_ENUM(Problem)
RENAME_ENUM(RecycleSubspace)
RENAME_ENUM(ConjType)
RENAME_ENUM(NO_CONJ)
RENAME_ENUM(CONJ)
RENAME_ENUM(MsgType)
RENAME_ENUM(Errors)
RENAME_ENUM(Warnings)
RENAME_ENUM(IterationDetails)
RENAME_ENUM(OrthoDetails)
RENAME_ENUM(FinalSummary)
RENAME_ENUM(TimingDetails)
RENAME_ENUM(StatusTestDetails)
RENAME_ENUM(Debug)
#undef RENAME_ENUM

#define BELOS_DEPRECATED
%include "Belos_Types.i"
