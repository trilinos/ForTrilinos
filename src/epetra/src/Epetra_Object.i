%module Epetra_Object

%{
#include "Epetra_DLLExportMacro.h"
#include "Epetra_Object.h"
%}

// ================= TODO =================
%ignore Epetra_Object::GetTracebackStream;                  // return std::ostream&
%ignore Epetra_Object::Print;                               // takes in std::ostream&
%ignore Epetra_Object::ReportError;                         // takes in std::string
%ignore operator<<(std::ostream& os, const Epetra_Object&); // not a class member
// ========================================
/* %ignore Epetra_Object::SetLabel;                            // const char* const parameter */
%ignore Epetra_Object::Label;                               // returns const char*

%include "Epetra_DLLExportMacro.h"
%include "Epetra_Object.h"
