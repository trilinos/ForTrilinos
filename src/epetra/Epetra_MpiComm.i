%module Epetra

%{
#include "Epetra_DLLExportMacro.h"
#include "Epetra_Comm.h"
%}

// Fortran's INTEGER(C_LONG) is the same as either INTEGER(C_INT) or INTEGER(C_LONG_LONG)
// This breaks generic interfaces due to ambiguity
%ignore Epetra_Comm::Broadcast(long *, int, int) const;
%ignore Epetra_Comm::GatherAll(long *, long *, int) const;
%ignore Epetra_Comm::SumAll(long *, long *, int) const;
%ignore Epetra_Comm::MaxAll(long *, long *, int) const;
%ignore Epetra_Comm::MinAll(long *, long *, int) const;
%ignore Epetra_Comm::ScanSum(long *, long *, int) const;
%ignore Epetra_Comm::Clone;

// Don't need these for now
%ignore PrintInfo;
%ignore CreateDistributor;
%ignore CreateDirectory;

%include "Epetra_DLLExportMacro.h"
%include "Epetra_Comm.h"
