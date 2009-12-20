module FEpetra_Map
  use ForTrilinos_enums!,only: FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Object_ID_t
  use FEpetra_Comm, only: epetra_comm
  use FEpetra_BlockMap ,only: epetra_BlockMap
  use iso_c_binding ,only: c_int
  use forepetra
  private                     ! Hide everything by default
  public :: epetra_map ! Expose type/constructors/methods
  public :: epetra_BlockMap
  implicit none

  type, extends(epetra_BlockMap) :: epetra_map_object
    private 
    type(FT_Epetra_Map_ID_t) :: map_id
  end type

  type epetra_map_cocoon       !"construction shell"
    private
    type(epetra_map_object), pointer :: object => null()
  end type 

  type epetra_map !"shell"
    private
    type(epetra_map_object), pointer :: object => null()
  contains
     procedure :: get_ID => EpetraMap_ID 
     procedure :: generalize 
     procedure :: force_finalization 
     procedure :: assign_shell
     procedure :: assign_cocoon
     generic :: assignment(=) => assign_cocoon,assign_shell
     final :: finalize
  end type

   interface epetra_map ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface
 
contains

  type(FT_Epetra_Map_ID_t) function EpetraMap_ID(this)
    class(epetra_map), intent(in) :: this 
    if (associated(this%object)) then
     EpetraMap_ID=this%object%map_id
    else
     stop 'EpetraMap_ID is unassociated cannot give back an ID'
    end if
  end function
 
  type(epetra_map_cocoon) function from_struct(id)
     type(FT_Epetra_Map_ID_t) ,intent(in) :: id
     allocate(from_struct%object)
     from_struct%object%map_id = id
  end function

  ! Original C++ prototype:
  ! Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Create ( int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(epetra_map_cocoon) function from_scratch(NGlobalElements,IndexBase,comm)
   ! use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: NGlobalElements
    integer(c_int) ,intent(in) :: IndexBase
    class(epetra_comm)         :: comm
    type(FT_Epetra_Comm_ID_t)  :: commID 
    allocate(from_scratch%object)
    commID=comm%get_ID()
    from_scratch%object%map_id = Epetra_Map_Create(NGlobalElements,IndexBase,commID)
  end function

  ! Original C++ prototype:
  ! Epetra_Map(const Epetra_Map& map);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );

  type(epetra_map_cocoon) function duplicate(original)
    type(epetra_map) ,intent(in) :: original
    allocate(duplicate%object)
    duplicate%object%map_id = Epetra_Map_Duplicate(original%object%map_id)
  end function

  subroutine assign_cocoon(lhs,rhs)
    class(epetra_map), intent(out) :: lhs
    type(epetra_map_cocoon), intent(in) :: rhs
    lhs%object => rhs%object
  end subroutine

  subroutine assign_shell(lhs,rhs)
    class(epetra_map), intent(out) :: lhs
    type(epetra_map), intent(in) :: rhs
    stop 'In epetra_map%assign_shell: no copying allowed'
  end subroutine

  subroutine finalize(this)
    type(epetra_map) :: this
    call Epetra_Map_Destroy( this%object%map_id ) 
  end subroutine

  subroutine force_finalization(this)
    class(epetra_map) ,intent(inout) :: this
    if (associated(this%object)) then
      call finalize(this) 
      deallocate(this%object)
    else
      print *,' finalization for epetra_map received object with unassociated selfID'
    end if
  end subroutine

  type(ForTrilinos_Object_ID_t) function generalize(this)
    class(epetra_map) ,intent(in) ,target :: this

    generalize = Epetra_Map_Abstract ( this%object%map_id)

   !Alternate implementation:
   !use iso_c_binding ,only : c_loc
   !interface
   !  type(ForTrilinos_Object_ID_t) function FEpetra_Map_Abstract(object_id_ptr) bind(C,name="for_linking_only")
   !    use iso_c_binding ,only : c_ptr
   !    import :: ForTrilinos_Object_ID_t
   !    type(c_ptr) ,value :: object_id_ptr
   !  end function
   !end interface
   !epetra_serialcomm_generalize = FEpetra_Map_Abstract( c_loc(this%object%map_id) )

  end function
 
end module 

