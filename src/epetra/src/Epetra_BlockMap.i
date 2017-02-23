%module Epetra_BlockMap

%{
#include "Epetra_DLLExportMacro.h"
#include "Epetra_BlockMap.h"
%}

// POSTPONE
%ignore Epetra_BlockMap::MyGlobalElementsPtr;
%ignore Epetra_BlockMap::MyGlobalElements;
%ignore Epetra_BlockMap::DataPtr;
%ignore Epetra_BlockMap::ElementSizeList;

%include "Epetra_DLLExportMacro.h"
%include "Epetra_BlockMap.h"
