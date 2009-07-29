module FEpetra_Comm
use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t
implicit none

type ,extends(universal) :: epetra_comm
  private
  type(FT_Epetra_Comm_ID_t) :: selfID
contains
  procedure :: clone
  procedure :: generalize => epetra_comm_generalize
  procedure :: final_subroutine => finalize
  final :: finalize
end type

interface epetra_comm_ ! We can remove trailing underscore on Fortran 2003-compliant compilers
  procedure constructor
end interface

contains
  
  type(epetra_comm) function constructor(comm)
    class(universal) ,intent(in) :: comm
    constructor%selfID = Epetra_Comm_Cast( comm%generalize() )
  end function

  ! Original C++ prototype:
  ! virtual Epetra_Comm * Clone() const = 0;
  ! CTrilinos prototype:
  ! CT_Epetra_Comm_ID_t Epetra_Comm_Clone ( CT_Epetra_Comm_ID_t selfID );

  type(epetra_comm) function clone(this)
    class(epetra_comm) ,intent(in) :: this
    clone%selfID = Epetra_Comm_Clone(this%selfID)
  end function

  ! Original C++ prototype:
  ! virtual ~Epetra_Comm();
  ! CTrilinos prototype:
  ! void Epetra_Comm_Destroy ( CT_Epetra_Comm_ID_t * selfID );

  subroutine finalize(this)
    type(epetra_comm) :: this
    call Epetra_Comm_Destroy ( this%selfID ) 
  end subroutine

end module 
