module Epetra_Map_Module
  use ,intrinsic :: iso_c_binding ,only : c_int
  use :: forepetraext
  implicit none
  private

  type ,bind(C) :: Epetra_Map 
    private
    integer(c_int) :: id
  end type 

  public :: Epetra_Map ,Create ,NumGlobalElements
contains

  function Create(numGlobalElements) 
    integer(c_int) :: numGlobalElements
    type(Epetra_Map) :: Create
    Create%id = FEpetra_Map_Create(numGlobalElements)
  end function create
  
  function NumGlobalElements(map) result(numGlobalElements_rtn)
    type(Epetra_Map) :: map
    integer(c_int) :: numGlobalElements_rtn
    numGlobalElements_rtn = FEpetra_Map_NumGlobalElements(map%id)
  end function NumGlobalElements
  
end module Epetra_Map_Module
