module FEpetra_Comm
use ForTrilinos_enums ,only : FT_Epetra_SerialComm_ID_t
implicit none

type ,extends(universal) :: epetra_serialcomm
  private
  type(FT_Epetra_Comm_ID_t) :: selfID
contains
  procedure :: clone
  procedure :: generalize => epetra_serialcomm_generalize
  procedure :: final_subroutine => finalize
  final :: finalize
end type

interface epetra_serialcomm ! constructors
  procedure from_scratch,duplicate
end interface

contains

  ! Original C++ prototype:
  ! Epetra_SerialComm();
  ! CTrilinos prototype:
  ! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  );
  
  type(epetra_serialcomm) function from_scratch()
    from_scratch%selfID = Epetra_SerialComm_Create()
  end function

  ! Original C++ prototype:
  ! Epetra_SerialComm(const Epetra_SerialComm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( CT_Epetra_SerialComm_ID_t CommID );

  type(epetra_serialcomm) function duplicate(original)
    type(epetra_serialcomm) ,intent(in) :: original
    duplicate%selfID = Epetra_SerialComm_Duplicate(original%selfID)
  end function

  ! Original C++ prototype:
  ! Epetra_Comm * Clone() const;
  ! CTrilinos prototype:
  ! CT_Epetra_Comm_ID_t Epetra_SerialComm_Clone ( CT_Epetra_SerialComm_ID_t selfID );
  
  type(epetra_comm) function clone(this)
    class(epetra_serialcomm) ,intent(in) :: this
    clone%selfID = Epetra_SerialComm_Clone(this%selfID)
  end function

  ! Original C++ prototype:
  ! virtual ~Epetra_SerialComm();
  ! CTrilinos prototype:
  ! void Epetra_SerialComm_Destroy ( CT_Epetra_SerialComm_ID_t * selfID );

  subroutine finalize(this)
    type(epetra_serialcomm) :: this
    call Epetra_SerialComm_Destroy ( this%selfID ) 
  end subroutine

end module 
