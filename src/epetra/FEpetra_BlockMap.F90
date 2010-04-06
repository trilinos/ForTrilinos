module FEpetra_BlockMap
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal
  use FEpetra_Comm  ,only: epetra_comm
  use iso_c_binding ,only: c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: epetra_BlockMap ! Expose type/constructors/methods

  type ,extends(universal)                 :: epetra_BlockMap !"shell"
    private
    type(FT_Epetra_BlockMap_ID_t) ,pointer :: BlockMap_id => null()
  contains
     !Developers only
     procedure         :: get_EpetraBlockMap_ID 
     procedure ,nopass :: alias_EpetraBlockMap_ID
     procedure         :: generalize 
     procedure         :: assign_to_epetra_BlockMap
     generic :: assignment(=) => assign_to_epetra_BlockMap
     !Local/Global ID accessor methods
     !Size and dimension acccessor functions
     procedure         :: NumGlobalElements
     procedure         :: NumMyElements
     procedure         :: MyGlobalElements
     procedure         :: ElementSize_Const
     procedure         :: ElementSize_LID
     generic :: ElementSize=>ElementSize_Const,ElementSize_LID
     !Miscellaneous boolean tests
     !Array accessor functions
     !Miscellaneous
     !Memory Management
     procedure         :: force_finalization 
     final :: finalize
  end type

   interface epetra_BlockMap ! constructors
     module procedure from_scratch,duplicate,from_struct,from_scratch_linear,from_scratch_arbitrary,from_scratch_variable
   end interface
 
contains
  type(FT_Epetra_BlockMap_ID_t) function from_struct(id)
     type(FT_Epetra_BlockMap_ID_t) ,intent(in) :: id
     from_struct = id
  end function
 
  ! Original C++ prototype:
  ! Epetra_BlockMap(int NumGlobalElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create ( int NumGlobalElements, int ElementSize, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_BlockMap_ID_t) function from_scratch(Num_GlobalElements,ElementSize,IndexBase,comm)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: ElementSize
    integer(c_int) ,intent(in) :: IndexBase
    class(epetra_comm)         :: comm
    from_scratch = Epetra_BlockMap_Create(Num_GlobalElements,ElementSize,IndexBase,comm%get_EpetraComm_ID())
  end function

! Original C++ prototype:
  ! Epetra_BlockMap(int NumGlobalElements, int NumMyElements, int ElementSize, int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Linear ( int NumGlobalElements, int NumMyElements, int
  ! ElementSize, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_BlockMap_ID_t) function from_scratch_linear(Num_GlobalElements,Num_MyElements,Element_Size,IndexBase,comm)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) :: Element_Size
    integer(c_int) ,intent(in) :: IndexBase
    class(epetra_comm)         :: comm
    from_scratch_linear = Epetra_BlockMap_Create_Linear(Num_GlobalElements,Num_MyElements,Element_Size,IndexBase,comm%get_EpetraComm_ID())
  end function

!Original C++ prototype:
  ! Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, int ElementSize, int IndexBae, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Arbitrary ( int NumGlobalElements, int NumMyElements,
  !const int * MyGlobalElements, int ElementSize, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_BlockMap_ID_t) function from_scratch_arbitrary(Num_GlobalElements,Num_MyElements,My_GlobalElements,Element_Size,IndexBase,comm)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements  
    integer(c_int) ,intent(in) :: Element_Size
    integer(c_int) ,intent(in) :: IndexBase
    class(epetra_comm)         :: comm
    from_scratch_arbitrary = Epetra_BlockMap_Create_Arbitrary(Num_GlobalElements,Num_MyElements,My_GlobalElements,Element_Size,IndexBase,comm%get_EpetraComm_ID())
  end function

! Original C++ prototype:
  ! Epetra_BlockMap(int NumGlobalElements, int NumMyElements, const int *MyGlobalElements, const int *ElementSizeList,
  !int IndexBase, const Epetra_Comm& Comm);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Create_Variable ( int NumGlobalElements, int NumMyElements,
  !const int * MyGlobalElements, const int * ElementSizeList, int IndexBase, CT_Epetra_Comm_ID_t CommID );

  type(FT_Epetra_BlockMap_ID_t) function from_scratch_variable(Num_GlobalElements,Num_MyElements,My_GlobalElements,Element_SizeList,IndexBase,comm)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_Map_ID_t
    integer(c_int) ,intent(in) :: Num_GlobalElements
    integer(c_int) ,intent(in) :: Num_MyElements
    integer(c_int) ,intent(in) ,dimension(:) :: My_GlobalElements  
    integer(c_int) ,intent(in) ,dimension(:) :: Element_SizeList    
    integer(c_int) ,intent(in) :: IndexBase
    class(epetra_comm)         :: comm
    from_scratch_variable = Epetra_BlockMap_Create_Variable(Num_GlobalElements,Num_MyElements,My_GlobalElements,Element_SizeList,IndexBase,comm%get_EpetraComm_ID())
  end function

  ! Original C++ prototype:
  ! Epetra_BlockMap(const Epetra_BlockMap& map);
  ! CTrilinos prototype:
  ! CT_Epetra_BlockMap_ID_t Epetra_BlockMap_Duplicate ( CT_Epetra_BlockMap_ID_t mapID );

  type(FT_Epetra_BlockMap_ID_t) function duplicate(original)
    type(epetra_BlockMap) ,intent(in) :: original
    duplicate = Epetra_BlockMap_Duplicate(original%BlockMap_id)
  end function

  type(FT_Epetra_BlockMap_ID_t) function get_EpetraBlockMap_ID(this)
    class(epetra_BlockMap) ,intent(in) :: this 
    if (associated(this%BlockMap_id)) then
     get_EpetraBlockMap_ID=this%BlockMap_id
    else
     stop 'get_EpetraBlockMap_ID: BlockMap_id is unassociated'
    end if
  end function
  
  type(FT_Epetra_BlockMap_ID_t) function alias_EpetraBlockMap_ID(generic_id)
    use ForTrilinos_table_man
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Epetra_BlockMap_ID
    use iso_c_binding     ,only: c_loc
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_BlockMap_ID))
    alias_EpetraBlockMap_ID=degeneralize_EpetraBlockMap(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(epetra_BlockMap) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%BlockMap_ID))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(epetra_BlockMap) ,intent(in) ,target :: this
   !generalize = Epetra_BlockMap_Generalize ( this%BlockMap_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_BlockMap_ID_t) function degeneralize_EpetraBlockMap(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_BlockMap_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                   ,value   :: generic_id
    type(FT_Epetra_BlockMap_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraBlockMap = local_ptr
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
   !degeneralize_EpetraBlockMap = Epetra_BlockMap_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  subroutine assign_to_epetra_BlockMap(lhs,rhs)
    class(epetra_BlockMap)        ,intent(inout) :: lhs
    type(FT_Epetra_BlockMap_ID_t) ,intent(in)    :: rhs
    allocate(lhs%BlockMap_id,source=rhs)
  end subroutine

  integer(c_int) function NumGlobalElements(this)
    class(epetra_BlockMap) ,intent(in) :: this
    NumGlobalElements=Epetra_BlockMap_NumGlobalElements(this%BlockMap_id)
  end function 

  integer(c_int) function NumMyElements(this)
    class(epetra_BlockMap) ,intent(in) :: this
    NumMyElements=Epetra_BlockMap_NumMyElements(this%BlockMap_id)
  end function 
 
  function MyGlobalElements(this) result(MyGlobalElementsList)
    class(epetra_BlockMap)     ,intent(in)    :: this
    integer(c_int),dimension(:),allocatable   :: MyGlobalElementsList
    integer(c_int)                            :: junk
    allocate(MyGlobalElementsList(this%NumMyElements()))
    junk=Epetra_BlockMap_MyGlobalElements_Fill(this%BlockMap_id,MyGlobalElementsList)
  end function 

  integer(c_int) function ElementSize_Const(this)
    class(epetra_BlockMap) ,intent(in) :: this
    ElementSize_Const=Epetra_BlockMap_ElementSize_Const(this%BlockMap_id)
  end function 

  integer(c_int) function ElementSize_LID(this,L_ID)
    class(epetra_BlockMap) ,intent(in) :: this
    integer(c_int)         ,intent(in) :: L_ID
    ElementSize_LID=Epetra_BlockMap_ElementSize(this%BlockMap_id,L_ID)
  end function 

  subroutine finalize(this)
    type(epetra_BlockMap) :: this
    print *,'finalize_BlockMap'
    call Epetra_BlockMap_Destroy( this%BlockMap_id ) 
    deallocate (this%BlockMap_id)
  end subroutine

  subroutine force_finalization(this)
    class(epetra_BlockMap) ,intent(inout) :: this
    if (associated(this%BlockMap_id)) then
      call finalize(this) 
    else
      print *,' finalization for epetra_BlockMap received  with unassociated object'
    end if
  end subroutine

end module 

