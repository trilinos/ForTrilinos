module FEpetra_Map
  use ForTrilinos_enums !,only: FT_Epetra_BlockMap_ID,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use FEpetra_Comm       ,only: epetra_comm
  use FEpetra_BlockMap   ,only: epetra_BlockMap
  use iso_c_binding      ,only: c_int
  use forepetra
  private                     ! Hide everything by default
  public :: epetra_map,epetra_BlockMap ! Expose type/constructors/methods
  implicit none

  type, extends(epetra_BlockMap)      :: epetra_map !"shell"
    private
    type(FT_Epetra_Map_ID_t) ,pointer :: map_id => null()
  contains
     !Developers only
     procedure         :: get_EpetraMap_ID 
     procedure ,nopass :: alias_EpetraMap_ID
     procedure         :: generalize 
     procedure         :: assign_to_epetra_Map
     procedure         :: epetraMap_assign_to_epetraMap  
     generic :: assignment(=) => assign_to_epetra_Map,epetraMap_assign_to_epetraMap
     !Memory Management
     procedure         :: force_finalization 
     final :: finalize
  end type

   interface epetra_map ! constructors
     module procedure from_scratch,duplicate,from_struct,from_scratch_linear,from_scratch_arbitrary
   end interface
 
contains

  type(FT_Epetra_Map_ID_t) function from_struct(id)
     type(FT_Epetra_Map_ID_t) ,intent(in) :: id
     from_struct = id
  end function

  ! Original C++ prototype:
  ! Epetra_Map(int NumGlobalElements, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Create ( int NumGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_Map_ID_t) function from_scratch(Num_GlobalElements,IndexBase,comm)
    use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: IndexBase
    class(epetra_comm)         :: comm
    type(FT_Epetra_Comm_ID_t) :: test   ! test line
    PRINT *,'start create map'          ! test line
    test=comm%get_EpetraComm_ID()       ! test line
    PRINT *,test%table,test%index,test%is_const ! test line
    from_scratch = Epetra_Map_Create(Num_GlobalElements,IndexBase,comm%get_EpetraComm_ID())
    PRINT *,'end create map'
  end function

! Original C++ prototype:
  ! Epetra_Map(int NumGlobalElements, int NumMyElements, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Create_Linear ( int NumGlobalElements, int NumMyElements, int IndexBase,
  ! CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_Map_ID_t) function from_scratch_linear(Num_GlobalElements,Num_MyElements,IndexBase,comm)
    use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) :: IndexBase
    class(epetra_comm)         :: comm
    from_scratch_linear = Epetra_Map_Create_Linear(Num_GlobalElements,Num_MyElements,IndexBase,comm%get_EpetraComm_ID())
  end function

! Original C++ prototype:
  ! Epetra_Map(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int IndexBase, const
  ! Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Create_Arbitrary ( int NumGlobalElements, int NumMyElements, const int
  !* MyGlobalElements, int IndexBase, CT_Epetra_Comm_ID_t CommID );

 type(FT_Epetra_Map_ID_t) function from_scratch_arbitrary(Num_GlobalElements,Num_MyElements,My_GlobalElements,IndexBase,comm)
    use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in)              :: Num_GlobalElements
    integer(c_int) ,intent(in)              :: Num_MyElements
    integer(c_int) ,intent(in) ,dimension(:),allocatable:: My_GlobalElements
    integer(c_int) ,intent(in)              :: IndexBase
    class(epetra_comm)                      :: comm
    from_scratch_arbitrary = Epetra_Map_Create_Arbitrary(Num_GlobalElements,Num_MyElements,My_GlobalElements,IndexBase,comm%get_EpetraComm_ID())
  end function

  ! Original C++ prototype:
  ! Epetra_Map(const Epetra_Map& map);
  ! CTrilinos prototype:
  ! CT_Epetra_Map_ID_t Epetra_Map_Duplicate ( CT_Epetra_Map_ID_t mapID );

  type(FT_Epetra_Map_ID_t) function duplicate(original)
    type(epetra_map) ,intent(in) :: original
    duplicate = Epetra_Map_Duplicate(original%map_id)
  end function

  type(FT_Epetra_Map_ID_t) function get_EpetraMap_ID(this)
    class(epetra_map) ,intent(in) :: this 
    if (associated(this%map_id)) then
     get_EpetraMap_ID=this%map_id
    else
     stop 'get_EpetraMap_ID: map_id is unassociated'
    end if
  end function
 
  type(FT_Epetra_Map_ID_t) function alias_EpetraMap_ID(generic_id)
    use iso_c_binding        ,only: c_loc
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t, FT_Epetra_Map_ID
    use ForTrilinos_table_man,only: CT_Alias
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Map_ID))
    alias_EpetraMap_ID=degeneralize_EpetraMap(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all 
   use iso_c_binding     ,only: c_loc
   class(epetra_map) ,intent(in) ,target :: this
   generalize =generalize_all(c_loc(this%map_id))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(epetra_map) ,intent(in) ,target :: this
   !generalize = Epetra_Map_Generalize ( this%map_id )
   ! ____ Use for CTrilinos function implementation ______
  
  end function
  
  type(FT_Epetra_Map_ID_t) function degeneralize_EpetraMap(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_Map_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)              ,value   :: generic_id
    type(FT_Epetra_Map_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraMap = local_ptr
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! type(ForTrilinos_Universal_ID_t) ,intent(in) : generic_id
   ! degeneralize_EpetraMap = Epetra_Map_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine assign_to_epetra_Map(lhs,rhs)
    class(epetra_map)        ,intent(inout):: lhs
    type(FT_Epetra_Map_ID_t) ,intent(in)   :: rhs
    allocate(lhs%map_id,source=rhs)
    lhs%epetra_BlockMap=epetra_BlockMap(lhs%alias_EpetraBlockMap_ID(lhs%generalize()))
  end subroutine
  
  subroutine epetraMap_assign_to_epetraMap(lhs,rhs)
    class(epetra_map) ,intent(inout):: lhs
    class(epetra_map) ,intent(in)   :: rhs
    call Epetra_Map_Assign(lhs%map_id,rhs%map_id)
    lhs = epetra_map(lhs%map_id)
  end subroutine

  subroutine finalize(this)
    type(epetra_map) :: this
    call Epetra_Map_Destroy( this%map_id ) 
    deallocate(this%map_id)
  end subroutine

  subroutine force_finalization(this)
    class(epetra_map) ,intent(inout) :: this
    if (associated(this%map_id)) then
      call finalize(this) 
    else
      print *,' finalization for epetra_map received object with unassociated map_id'
    end if
  end subroutine

end module 

