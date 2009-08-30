module FEpetra_SerialComm
  use ForTrilinos_enums ,only : FT_Epetra_SerialComm_ID_t,ForTrilinos_Object_ID_t
  use ForTrilinos_universal ,only : universal
  use forepetra
  private                     ! Hide everything by default
  public :: epetra_serialcomm ! Expose type/constructor/methods
  implicit none

  type ,extends(universal) :: epetra_serialcomm
    private
    type(FT_Epetra_SerialComm_ID_t) :: selfID
  contains
    procedure :: clone
    procedure :: generalize => epetra_serialcomm_generalize
    procedure :: force_finalization => call_finalize
    final :: finalize
  end type

  interface epetra_serialcomm ! constructors
    procedure from_scratch,duplicate
  end interface

contains

  ! CTrilinos prototype:
  ! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Cast ( CTrilinos_Object_ID_t id );


  ! CTrilinos prototype:
  ! CTrilinos_Object_ID_t Epetra_SerialComm_Abstract ( CT_Epetra_SerialComm_ID_t id );

  type(ForTrilinos_Object_ID_t) function epetra_serialcomm_generalize(this)
    use iso_c_binding ,only : c_loc
    class(epetra_serialcomm) ,intent(in) ,target :: this

   !interface
   !  type(ForTrilinos_Object_ID_t) function FEpetra_SerialComm_Abstract(object_id_ptr) bind(C,name="for_linking_only")
   !    use iso_c_binding ,only : c_ptr
   !    import :: ForTrilinos_Object_ID_t
   !    type(c_ptr) ,value :: object_id_ptr
   !  end function
   !end interface

    epetra_serialcomm_generalize = Epetra_SerialComm_Abstract ( this%selfID )
   !epetra_serialcomm_generalize = FEpetra_SerialComm_Abstract( c_loc(this%selfID) )
  end function

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
    use FEpetra_Comm ,only : epetra_comm
    class(epetra_serialcomm) ,intent(in) :: this
    clone = Epetra_SerialComm_Clone(this%selfID)
  end function

  ! Original C++ prototype:
  ! virtual ~Epetra_SerialComm();
  ! CTrilinos prototype:
  ! void Epetra_SerialComm_Destroy ( CT_Epetra_SerialComm_ID_t * selfID );

  subroutine finalize(this)
    type(epetra_serialcomm) :: this
    call Epetra_SerialComm_Destroy ( this%selfID ) 
  end subroutine

  subroutine call_finalize(this)
    class(epetra_serialcomm) ,intent(inout) :: this
    call finalize(this) 
  end subroutine

end module 
