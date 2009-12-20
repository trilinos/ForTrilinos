module FEpetra_SerialComm
  use ForTrilinos_enums !,only : FT_Epetra_SerialComm_ID_t,ForTrilinos_Object_ID_t
  use FEpetra_Comm ,only : epetra_comm
  use forepetra
  private                     ! Hide everything by default
  public :: epetra_serialcomm ! Expose type/constructors/methods
  public :: epetra_comm ! Expose type/constructors/methods
  implicit none

  type,extends(epetra_comm) :: epetra_serialcomm_object
    private
    type(FT_Epetra_SerialComm_ID_t) :: serialcomm_id 
  end type

  type epetra_serialcomm_cocoon !"construction shell"
    private
    type(epetra_serialcomm_object), pointer :: object => null()
  end type 

  type epetra_serialcomm        !"shell"
    private
    type(epetra_serialcomm_object), pointer :: object => null()
  contains
     procedure :: get_ID => EpetraSerialComm_ID 
     procedure :: generalize 
     procedure :: force_finalization 
!     procedure :: clone 
     procedure :: assign_shell
     procedure :: assign_cocoon
     generic :: assignment(=) => assign_cocoon,assign_shell
     final :: finalize
  end type

   interface epetra_serialcomm ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface

contains

 type(FT_Epetra_SerialComm_ID_t) function EpetraSerialComm_ID(this)
   class(epetra_serialcomm), intent(in) :: this 
       if (associated(this%object)) then 
         EpetraSerialComm_ID=this%object%serialcomm_id
       else
         stop 'EpetraSerialComm_object is unassociated cannot give back a ID'
       end if
  end function

  type(epetra_serialcomm_cocoon) function from_struct(id)
     type(FT_Epetra_SerialComm_ID_t) ,intent(in) :: id
     allocate(from_struct%object)
     from_struct%object%serialcomm_id = id
  end function

  ! Original C++ prototype:
  ! Epetra_SerialComm();
  ! CTrilinos prototype:
  ! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Create (  );
  
  type(epetra_serialcomm_cocoon) function from_scratch()
    allocate(from_scratch%object)
    from_scratch%object%serialcomm_id = Epetra_SerialComm_Create()
  end function

  ! Original C++ prototype:
  ! Epetra_SerialComm(const Epetra_SerialComm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_SerialComm_ID_t Epetra_SerialComm_Duplicate ( CT_Epetra_SerialComm_ID_t CommID );

  type(epetra_serialcomm_cocoon) function duplicate(original)
    type(epetra_serialcomm) ,intent(in) :: original
    allocate(duplicate%object)
    duplicate%object%serialcomm_id = Epetra_SerialComm_Duplicate(original%object%serialcomm_id)
  end function

!  function clone(this)
!    class(epetra_serialcomm) ,intent(in) :: this
!    class(epetra_comm) ,allocatable :: clone
!    allocate(epetra_serialcomm :: clone) 
!    select type(clone)
!      class is (epetra_comm)
!       !!clone%comm_id = Epetra_SerialComm_Clone(this%object%serialcomm_id)
!       !! stop 'epetra_serialcomm%clone(): determine return type now that epetra_comm is abstract'
!      class default
!        stop 'epetra_serialcomm%clone(): class not supported'
!    end select
!  end function

  subroutine assign_cocoon(lhs,rhs)
    class(epetra_serialcomm), intent(out) :: lhs
    type(epetra_serialcomm_cocoon), intent(in) :: rhs
    lhs%object => rhs%object
  end subroutine

  subroutine assign_shell(lhs,rhs)
    class(epetra_serialcomm), intent(out) :: lhs
    type(epetra_serialcomm), intent(in) :: rhs
    stop 'In epetra_serialcomm%assign_shell: no copying allowed'
  end subroutine

  subroutine finalize(this)
    type(epetra_serialcomm) :: this
    call Epetra_SerialComm_Destroy( this%object%serialcomm_id ) 
  end subroutine

  subroutine force_finalization(this)
    class(epetra_serialcomm) ,intent(inout) :: this
    if (associated(this%object)) then
      call finalize(this) 
      deallocate(this%object)
    else
      print *,' finalization for epetra_serialcomm received object with unassociated serialcomm_id'
    end if
  end subroutine

  type(ForTrilinos_Object_ID_t) function generalize(this)
    class(epetra_serialcomm) ,intent(in) ,target :: this
    
    generalize = Epetra_SerialComm_Abstract ( this%object%serialcomm_id )

   !Alternate implementation:
   !use iso_c_binding ,only : c_loc
   !interface
   !  type(ForTrilinos_Object_ID_t) function FEpetra_SerialComm_Abstract(object_id_ptr) bind(C,name="for_linking_only")
   !    use iso_c_binding ,only : c_ptr
   !    import :: ForTrilinos_Object_ID_t
   !    type(c_ptr) ,value :: object_id_ptr
   !  end function
   !end interface
   !epetra_serialcomm_generalize = FEpetra_SerialComm_Abstract( c_loc(this%serialcomm_id) )

  end function
 
end module 

