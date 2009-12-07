module FEpetra_SerialComm
  use ForTrilinos_enums ,only : FT_Epetra_SerialComm_ID_t,ForTrilinos_Object_ID_t
  use FEPetra_Comm ,only : Epetra_Comm
  use forepetra
  private                     ! Hide everything by default
  public :: epetra_serialcomm ! Expose type/constructors/methods
  implicit none

  type ,extends(Epetra_Comm) :: epetra_serialcomm
    private
    type(FT_Epetra_SerialComm_ID_t) :: selfID
  contains
     procedure :: generalize 
     procedure :: force_finalization 
     procedure :: clone 
     final :: finalize
  end type

   interface epetra_serialcomm ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface
 
contains

  type(epetra_serialcomm) function from_struct(id)
     type(FT_Epetra_SerialComm_ID_t) ,intent(in) :: id
     from_struct%selfID = id
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

  function clone(this)
    class(epetra_serialcomm) ,intent(in) :: this
    class(epetra_comm) ,allocatable :: clone
    allocate(epetra_serialcomm :: clone) 
    select type(clone)
      class is (epetra_serialcomm)
       !clone%selfID = Epetra_SerialComm_Clone(this%selfID)
        stop 'epetra_serialcomm%clone(): determine return type now that epetra_comm is abstract'
      class default
        stop 'epetra_serialcomm%clone(): class not supported'
    end select
  end function

  subroutine finalize(this)
    type(epetra_serialcomm) :: this
    call Epetra_SerialComm_Destroy( this%selfID ) 
  end subroutine

  subroutine force_finalization(this)
    class(epetra_serialcomm) ,intent(inout) :: this
    call finalize(this) 
  end subroutine

  type(ForTrilinos_Object_ID_t) function generalize(this)
    class(epetra_serialcomm) ,intent(in) ,target :: this

    generalize = Epetra_SerialComm_Abstract ( this%selfID )

   !Alternate implementation:
   !use iso_c_binding ,only : c_loc
   !interface
   !  type(ForTrilinos_Object_ID_t) function FEpetra_SerialComm_Abstract(object_id_ptr) bind(C,name="for_linking_only")
   !    use iso_c_binding ,only : c_ptr
   !    import :: ForTrilinos_Object_ID_t
   !    type(c_ptr) ,value :: object_id_ptr
   !  end function
   !end interface
   !epetra_serialcomm_generalize = FEpetra_SerialComm_Abstract( c_loc(this%selfID) )

  end function
 
end module 

