module FEpetra_MpiComm
use ForTrilinos_enums ,only : FT_Epetra_MpiComm_ID_t
implicit none

type ,extends(universal) :: epetra_mpicomm
  private
  type(FT_Epetra_Comm_ID_t) :: selfID
contains
  procedure :: clone
  procedure :: generalize => epetra_mpicomm_generalize
  procedure :: final_subroutine => finalize
  final :: finalize
end type

interface epetra_mpicomm ! constructors
  procedure from_scratch,duplicate
end interface

contains

  ! Original C++ prototype:
  ! Epetra_MpiComm();
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create (  );
  
  type(epetra_mpicomm) function from_scratch()
    from_scratch%selfID = Epetra_MpiComm_Create()
  end function

  ! Original C++ prototype:
  ! Epetra_MpiComm(const Epetra_MpiComm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( CT_Epetra_MpiComm_ID_t CommID );

  type(epetra_mpicomm) function duplicate(original)
    type(epetra_mpicomm) ,intent(in) :: original
    duplicate%selfID = Epetra_MpiComm_Duplicate(original%selfID)
  end function

  ! Original C++ prototype:
  ! Epetra_Comm * Clone() const;
  ! CTrilinos prototype:
  ! CT_Epetra_Comm_ID_t Epetra_MpiComm_Clone ( CT_Epetra_MpiComm_ID_t selfID );
  
  type(epetra_comm) function clone(this)
    class(epetra_mpicomm) ,intent(in) :: this
    clone%selfID = Epetra_MpiComm_Clone(this%selfID)
  end function

  ! Original C++ prototype:
  ! virtual ~Epetra_MpiComm();
  ! CTrilinos prototype:
  ! void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID );

  subroutine finalize(this)
    type(epetra_mpicomm) :: this
    call Epetra_MpiComm_Destroy ( this%selfID ) 
  end subroutine

end module 
