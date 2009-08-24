module FEpetra_MpiComm
#ifdef HAVE_MPI
  use ForTrilinos_enums ,only : FT_Epetra_MpiComm_ID_t,ForTrilinos_Object_ID_t
  use ForTrilinos_universal ,only : universal
  use forepetra
  implicit none
  
  type ,extends(universal) :: epetra_mpicomm
    private
    type(FT_Epetra_MpiComm_ID_t) :: selfID
  contains
   !procedure :: clone
    procedure :: generalize => epetra_mpicomm_generalize
    procedure :: force_finalization => call_finalize
    final :: finalize
  end type
  
 !interface epetra_mpicomm ! constructors
 !  procedure from_scratch,duplicate
 !end interface
  
contains

  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Cast ( CTrilinos_Object_ID_t id );


  ! CTrilinos prototype:
  ! CTrilinos_Object_ID_t Epetra_MpiComm_Abstract ( CT_Epetra_MpiComm_ID_t id );

  type(ForTrilinos_Object_ID_t) function epetra_mpicomm_generalize(this)
    use iso_c_binding ,only : c_loc
    use forepetra ,only : Epetra_MpiComm_Abstract
    class(epetra_mpicomm) ,intent(in) ,target :: this

    epetra_mpicomm_generalize = Epetra_MpiComm_Abstract( this%selfID )
                                
   !interface
   !  type(ForTrilinos_Object_ID_t) function FEpetra_MpiComm_Abstract(object_id_ptr) bind(C,name="for_linking_only")
   !    use iso_c_binding ,only : c_ptr
   !    import :: ForTrilinos_Object_ID_t
   !    type(c_ptr) ,value :: object_id_ptr
   !  end function
   !end interface

   !epetra_mpicomm_generalize = FEpetra_Comm_Abstract( c_loc(this%selfID) )
  end function


  ! Original C++ prototype:
  ! Epetra_MpiComm();
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Create (  );
  
 !type(epetra_mpicomm) function from_scratch()
 !  from_scratch%selfID = Epetra_MpiComm_Create()
 !end function

  ! Original C++ prototype:
  ! Epetra_MpiComm(const Epetra_MpiComm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_MpiComm_ID_t Epetra_MpiComm_Duplicate ( CT_Epetra_MpiComm_ID_t CommID );

 !type(epetra_mpicomm) function duplicate(original)
 !  type(epetra_mpicomm) ,intent(in) :: original
 !  duplicate%selfID = Epetra_MpiComm_Duplicate(original%selfID)
 !end function

  ! Original C++ prototype:
  ! Epetra_Comm * Clone() const;
  ! CTrilinos prototype:
  ! CT_Epetra_Comm_ID_t Epetra_MpiComm_Clone ( CT_Epetra_MpiComm_ID_t selfID );
  
 !type(epetra_mpicomm) function clone(this)
 !  class(epetra_mpicomm) ,intent(in) :: this
 !  clone%selfID = Epetra_MpiComm_Clone(this%selfID)
 !end function

  ! Original C++ prototype:
  ! virtual ~Epetra_MpiComm();
  ! CTrilinos prototype:
  ! void Epetra_MpiComm_Destroy ( CT_Epetra_MpiComm_ID_t * selfID );

  subroutine finalize(this)
    type(epetra_mpicomm) :: this
    call Epetra_MpiComm_Destroy ( this%selfID ) 
  end subroutine

  subroutine call_finalize(this)
    class(epetra_mpicomm) ,intent(inout) :: this
    call finalize(this) 
  end subroutine
#endif
end module 
