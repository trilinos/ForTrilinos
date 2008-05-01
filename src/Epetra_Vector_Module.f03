module Epetra_Vector_Module

  use ,intrinsic :: iso_c_binding     ,only : c_int
  use            :: epetra_map_module ,only : Epetra_Map, GetMapID
  use            :: forepetraext
  implicit none
  private

  type ,bind(C) ,public :: EPetra_Vector 
    private
    integer(c_int) :: id
  end type 

  interface Create
    module procedure Create_Epetra_Vector
  end interface 

  public :: Create
contains

  function Create_Epetra_Vector(map) 
    type(Epetra_Map) :: map
    type(Epetra_Vector) :: Create_Epetra_Vector
    Create_Epetra_Vector%id = FEpetra_Vector_Create(GetMapID(map))
  end function

end module Epetra_Vector_Module
