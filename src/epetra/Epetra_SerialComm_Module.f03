module Epetra_SerialComm_Module
  use ,intrinsic :: iso_c_binding ,only : c_int
  use            :: forepetra
  implicit none
  private
  public :: Epetra_SerialComm
  type ,public :: Epetra_SerialComm 
    private
    integer(c_int) :: id
  contains
    procedure :: Create
    procedure :: NumGlobalElements
    procedure :: GetSerialCommID
    procedure :: Destroy
   !final     :: Destroy
  end type 
  public :: Epetra_SerialComm_
 !interface Epetra_SerialComm
 !  procedure constructor
 !end interface 
contains
  ! Constructor
  type(Epetra_SerialComm) Epetra_SerialComm_() 
    Epetra_SerialComm_%id = FEpetra_SerialComm_Create()
  end subroutine
  
  ! Destructor 
  subroutine Destroy(serial_comm) 
    type(Epetra_SerialComm) ,intent(inout) :: serial_comm
    call FEpetra_SerialComm_Destroy(serial_comm%id)
  end subroutine 

  ! Accessors
  integer(c_int) function GetSerialCommID(serial_comm) 
    type(Epetra_SerialComm) ,intent(in) :: serial_comm
    GetSerialCommID = serial_comm%id
  end function 
end module Epetra_SerialComm_Module
