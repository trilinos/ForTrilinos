module Epetra_MpiComm_Module
  use ,intrinsic :: iso_c_binding ,only : c_int
  use            :: forepetra
  implicit none
  private
  public :: Epetra_MpiComm
  type ,public :: Epetra_MpiComm 
    private
    integer(c_int) :: id
  contains
    procedure :: Create
    procedure :: NumGlobalElements
    procedure :: GetMpiCommID
    procedure :: Destroy
   !final     :: Destroy
  end type 
  public :: Epetra_MpiComm_
 !interface Epetra_MpiComm
 !  procedure constructor
 !end interface 
contains
  ! Constructor
  type(Epetra_MpiComm) Epetra_MpiComm_() 
    Epetra_MpiComm_%id = FEpetra_MpiComm_Create()
  end subroutine
  
  ! Destructor 
  subroutine Destroy(mpi_comm) 
    type(Epetra_MpiComm) ,intent(inout) :: mpi_comm
    call FEpetra_MpiComm_Destroy(mpi_comm%id)
  end subroutine 

  ! Accessors
  integer(c_int) function GetMpiCommID(mpi_comm) 
    type(Epetra_MpiComm) ,intent(in) :: mpi_comm
    GetMpiCommID = mpi_comm%id
  end function 
end module Epetra_MpiComm_Module
