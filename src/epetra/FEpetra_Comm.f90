module FEpetra_Comm
use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t ,ForTrilinos_Object_ID_t
use ForTrilinos_universal ,only : universal
use forepetra
implicit none

type ,extends(universal) :: epetra_comm
  private
  type(FT_Epetra_Comm_ID_t) :: selfID
contains
  procedure :: clone
  procedure :: generalize => epetra_comm_generalize
  procedure :: invoke_final_subroutine => call_final
  final :: finalize
end type

interface epetra_comm
  procedure constructor
end interface

contains
  subroutine call_final(this)
    class(epetra_comm) ,intent(inout) :: this
    call finalize(this)
  end subroutine

  ! CTrilinos_Object_ID_t Epetra_Comm_Abstract (CT_Epetra_Comm_ID_t id );

  type(ForTrilinos_Object_ID_t) function epetra_comm_generalize(this) 
    class(epetra_comm) ,intent(in) :: this
    epetra_comm_generalize = Epetra_Comm_Abstract ( this%selfID )
  end function
  
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
