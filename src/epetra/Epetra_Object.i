%module Epetra_Object

%{
#include "Epetra_DLLExportMacro.h"
#include "Epetra_Object.h"
%}

// POSTPONE
%ignore Epetra_Object::Print;
%ignore Epetra_Object::ReportError;
%ignore Epetra_Object::GetTracebackStream;
%ignore operator<<(std::ostream& os, const Epetra_Object&);
%ignore Epetra_Object::SetLabel;
%ignore Epetra_Object::Label;
%ignore Epetra_Object::Epetra_Object; // disable constructors

%include "Epetra_DLLExportMacro.h"
%include "Epetra_Object.h"
