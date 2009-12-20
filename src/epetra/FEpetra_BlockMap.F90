module FEpetra_BlockMap
  use ForTrilinos_enums!,only: FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Object_ID_t
  use FEpetra_Comm, only: epetra_comm
  use ForTrilinos_universal
  use iso_c_binding ,only: c_int
  use forepetra
  private                     ! Hide everything by default
  public :: epetra_BlockMap ! Expose type/constructors/methods
  implicit none

  !type, extends(universal) :: epetra_BlockMap_object
  type epetra_BlockMap_object
    private 
    type(FT_Epetra_BlockMap_ID_t) :: BlockMap_id
  end type

  type epetra_BlockMap_cocoon       !"construction shell"
    private
    type(epetra_BlockMap_object), pointer :: object => null()
  end type 

  type epetra_BlockMap !"shell"
    private
    type(epetra_BlockMap_object), pointer :: object => null()
  contains
     procedure :: NumGlobalElements
     procedure :: get_ID => EpetraBlockMap_ID 
     procedure :: generalize 
     procedure :: force_finalization 
     procedure :: assign_shell
     procedure :: assign_cocoon
     generic :: assignment(=) => assign_cocoon,assign_shell
     final :: finalize
  end type

   interface epetra_BlockMap ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface
 
contains
  
  integer(c_int) function NumGlobalElements(this)
    class(epetra_BlockMap), intent(in) :: this
    NumGlobalElements=Epetra_BlockMap_NumGlobalElements(this%object%BlockMap_id)
  end function 

  type(FT_Epetra_BlockMap_ID_t) function EpetraBlockMap_ID(this)
    class(epetra_BlockMap), intent(in) :: this 
    if (associated(this%object)) then
     EpetraBlockMap_ID=this%object%BlockMap_id
    else
     stop 'EpetraBlockMap_ID is unassociated cannot give back an ID'
    end if
  end function
 
  type(epetra_BlockMap_cocoon) function from_struct(id)
     type(FT_Epetra_BlockMap_ID_t) ,intent(in) :: id
     allocate(from_struct%object)
     from_struct%object%BlockMap_id = id
  end function
 
  ! Original C++ prototype:
  ! Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( int NumGlobalElements, int ElementSize, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(epetra_BlockMap_cocoon) function from_scratch(NumGlobalElements,ElementSize,IndexBase,comm)
   ! use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: NumGlobalElements
    integer(c_int) ,intent(in) :: ElementSize
    integer(c_int) ,intent(in) :: IndexBase
    class(epetra_comm)         :: comm
    type(FT_Epetra_Comm_ID_t)  :: commID 
    allocate(from_scratch%object)
    commID=comm%get_ID()
    from_scratch%object%BlockMap_id = Epetra_BlockMap_Create(NumGlobalElements,ElementSize,IndexBase,commID)
  end function

  ! Original C++ prototype:
  ! Epetra_BlockMap(const Epetra_BlockMap& map);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( CT_Epetra_BlockMap_ID_t mapID );

  type(epetra_BlockMap_cocoon) function duplicate(original)
    type(epetra_BlockMap) ,intent(in) :: original
    allocate(duplicate%object)
    duplicate%object%BlockMap_id = Epetra_BlockMap_Duplicate(original%object%BlockMap_id)
  end function

  subroutine assign_cocoon(lhs,rhs)
    class(epetra_BlockMap), intent(out) :: lhs
    type(epetra_BlockMap_cocoon), intent(in) :: rhs
    lhs%object => rhs%object
  end subroutine

  subroutine assign_shell(lhs,rhs)
    class(epetra_BlockMap), intent(out) :: lhs
    type(epetra_BlockMap), intent(in) :: rhs
    stop 'In epetra_BlockMap%assign_shell: no copying allowed'
  end subroutine

  subroutine finalize(this)
    type(epetra_BlockMap) :: this
    call Epetra_BlockMap_Destroy( this%object%BlockMap_id ) 
  end subroutine

  subroutine force_finalization(this)
    class(epetra_BlockMap) ,intent(inout) :: this
    if (associated(this%object)) then
      call finalize(this) 
      deallocate(this%object)
    else
      print *,' finalization for epetra_BlockMap received object with unassociated object'
    end if
  end subroutine

  type(ForTrilinos_Object_ID_t) function generalize(this)
    class(epetra_BlockMap) ,intent(in) ,target :: this

    generalize = Epetra_BlockMap_Abstract ( this%object%BlockMap_id)

   !Alternate implementation:
   !use iso_c_binding ,only : c_loc
   !interface
   !  type(ForTrilinos_Object_ID_t) function FEpetra_BlockMap_Abstract(object_id_ptr) bind(C,name="for_linking_only")
   !    use iso_c_binding ,only : c_ptr
   !    import :: ForTrilinos_Object_ID_t
   !    type(c_ptr) ,value :: object_id_ptr
   !  end function
   !end interface
   !epetra_BlockMap_generalize = FEpetra_BlockMap_Abstract( c_loc(this%object%BlockMap_id) )

  end function
 
end module 

