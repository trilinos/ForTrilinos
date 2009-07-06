module Epetra_Map_Module
  use ,intrinsic :: iso_c_binding ,only : c_int
  use :: forepetra
  implicit none
  private

! type ,public :: Epetra_Map 
!   private
!   integer(c_int) :: id
! contains
!   procedure :: Create
!   procedure :: NumGlobalElements
!   procedure :: GetMapID
!   procedure :: Destroy
! end type 

contains

! subroutine Create(map,numGlobalElements) 
!   type(Epetra_Map) ,intent(out) :: map
!   integer(c_int)   ,intent(in)  :: numGlobalElements
!   map%id = FEpetra_Map_Create(numGlobalElements)
! end subroutine
  
! function GetMapID(map) 
!   type(Epetra_Map) ,intent(in) :: map
!   integer(c_int)               :: GetMapID
!   GetMapID = map%id
! end function 
  
! function NumGlobalElements(map)
!   type(Epetra_Map) ,intent(in) :: map
!   integer(c_int)               :: NumGlobalElements
!   NumGlobalElements = FEpetra_Map_NumGlobalElements(map%id)
! end function 

! subroutine Destroy(map) 
!   type(Epetra_Map) ,intent(inout) :: map
!   call FEpetra_Map_Destroy(map%id)
! end subroutine 
end module Epetra_Map_Module
