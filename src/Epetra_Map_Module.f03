module Epetra_Map_Module
  use ,intrinsic :: iso_c_binding ,only : c_int
  use :: forepetraext
  implicit none
  private

  type ,bind(C) ,public :: Epetra_Map 
    private
    integer(c_int) :: id
  end type 

  interface Create
    module procedure Create_Epetra_Map
  end interface

  public :: Create ,NumGlobalElements ,GetMapID
contains

  function Create_Epetra_Map(numGlobalElements) 
    integer(c_int) :: numGlobalElements
    type(Epetra_Map) :: Create_Epetra_Map
    Create_Epetra_Map%id = Epetra_Map_Create(numGlobalElements)
  end function 
  
  function GetMapID(map) 
    type(Epetra_Map) :: map
    integer(c_int)   :: GetMapID
    GetMapID = map%id
  end function 
  
  function NumGlobalElements(map) result(numGlobalElements_rtn)
    type(Epetra_Map) :: map
    integer(c_int) :: numGlobalElements_rtn
    numGlobalElements_rtn = Epetra_Map_NumGlobalElements(map%id)
  end function 
  
end module Epetra_Map_Module
