%module Epetra_SerialComm

%{
#include "Epetra_DLLExportMacro.h"
#include "Epetra_SerialComm.h"
%}

%include "Epetra_Object.i"
%include "Epetra_Comm.i"

// Fortran's INTEGER(C_LONG) is the same as either INTEGER(C_INT) or INTEGER(C_LONG_LONG)
// This breaks generic interfaces due to ambiguity
%ignore Epetra_SerialComm::Broadcast(long *, int, int) const;
%ignore Epetra_SerialComm::GatherAll(long *, long *, int) const;
%ignore Epetra_SerialComm::SumAll(long *, long *, int) const;
%ignore Epetra_SerialComm::MaxAll(long *, long *, int) const;
%ignore Epetra_SerialComm::MinAll(long *, long *, int) const;
%ignore Epetra_SerialComm::ScanSum(long *, long *, int) const;

// POSTPONE
%ignore Epetra_SerialComm::Clone;
%ignore Epetra_SerialComm::CreateDirectory;
%ignore Epetra_SerialComm::CreateDistributor;
%ignore Epetra_SerialComm::DataPtr;
%ignore Epetra_SerialComm::ReferenceCount;
%ignore Epetra_SerialComm::Print;
%ignore Epetra_SerialComm::PrintInfo;

%include "Epetra_DLLExportMacro.h"
%include "Epetra_SerialComm.h"
