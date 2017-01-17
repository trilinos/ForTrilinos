%module Epetra0

%{
#include "Epetra_DLLExportMacro.h"
%include "Epetra_Object.h"
%}

%ignore Epetra_Object::TracebackMode;
%ignore Epetra_Object::SetTracebackMode;
%ignore Epetra_Object::GetTracebackMode;
%ignore Epetra_Object::GetTracebackStream;


%include "Epetra_DLLExportMacro.h"
%include "Epetra_Object.h"
