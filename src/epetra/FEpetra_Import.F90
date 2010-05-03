module FEpetra_Import
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_Import_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal
  use FEpetra_Comm  ,only: Epetra_Comm
  use FEpetra_BlockMap ,only: Epetra_BlockMap
  use iso_c_binding ,only: c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: Epetra_Import ! Expose type/constructors/methods

  type ,extends(universal)                 :: Epetra_Import !"shell"
    private
    type(FT_Epetra_Import_ID_t) ,pointer :: Import_id => null()
  contains
     !Developers only
     procedure         :: get_EpetraImport_ID 
     procedure ,nopass :: alias_EpetraImport_ID
     procedure         :: generalize 
     procedure         :: assign_to_Epetra_Import
     generic :: assignment(=) => assign_to_Epetra_Import
     ! Public member functions
     procedure        :: NumSameIDs
     procedure        :: NumPermuteIDs
     !procedure        :: PermuteFromLIDs
     !procedure        :: PermuteToLIDs
     procedure        :: NumRemoteIDs
     !procedure       :: RemoteLIDs
     procedure        :: NumExportIDs
     !procedure       :: ExportLIDs
     !procedure       :: Export PIDs
     procedure        :: NumSend
     procedure        :: NumRecv
     procedure        :: SourceMap
     procedure        :: TargetMap
     !Memory Management
     procedure         :: force_finalization 
     final :: finalize
  end type

   interface Epetra_Import ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface
 
contains
  type(FT_Epetra_Import_ID_t) function from_struct(id)
     type(FT_Epetra_Import_ID_t) ,intent(in) :: id
     from_struct = id
  end function
 
  ! Original C++ prototype:
  ! Epetra_Import( const Epetra_BlockMap & TargetMap, const Epetra_BlockMap & SourceMap );
  ! CTrilinos prototype:
  ! CT_Epetra_Import_ID_t Epetra_Import_Create ( CT_Epetra_BlockMap_ID_t TargetMapID,
  ! CT_Epetra_BlockMap_ID_t SourceMapID );


  type(FT_Epetra_Import_ID_t) function from_scratch(TargetMap,SourceMap)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t
    type(Epetra_BlockMap), intent(in) :: TargetMap
    type(Epetra_BlockMap), intent(in) :: SourceMap
    from_scratch = Epetra_Import_Create(TargetMap%get_EpetraBlockMap_ID(),SourceMap%get_EpetraBlockMap_ID())
  end function

  ! Original C++ prototype:
  ! Epetra_Import(const Epetra_Import& Importer);
  ! CTrilinos prototype:
  ! CT_Epetra_Import_ID_t Epetra_Import_Duplicate ( CT_Epetra_Import_ID_t ImporterID );

  type(FT_Epetra_Import_ID_t) function duplicate(original)
    type(Epetra_Import) ,intent(in) :: original
    duplicate = Epetra_Import_Duplicate(original%Import_id)
  end function

  type(FT_Epetra_Import_ID_t) function get_EpetraImport_ID(this)
    class(Epetra_Import) ,intent(in) :: this 
    if (associated(this%Import_id)) then
     get_EpetraImport_ID=this%Import_id
    else
     stop 'get_EpetraImport_ID: Import_id is unassociated'
    end if
  end function
  
  type(FT_Epetra_Import_ID_t) function alias_EpetraImport_ID(generic_id)
    use ForTrilinos_table_man
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Epetra_Import_ID
    use iso_c_binding     ,only: c_loc
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Import_ID))
    alias_EpetraImport_ID=degeneralize_EpetraImport(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(Epetra_Import) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%Import_ID))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(Epetra_Import) ,intent(in) ,target :: this
   !generalize = Epetra_Import_Generalize ( this%Import_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_Import_ID_t) function degeneralize_EpetraImport(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_Import_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                   ,value   :: generic_id
    type(FT_Epetra_Import_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraImport = local_ptr
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
   !degeneralize_EpetraImport = Epetra_Import_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  subroutine assign_to_Epetra_Import(lhs,rhs)
    class(Epetra_Import)        ,intent(inout) :: lhs
    type(FT_Epetra_Import_ID_t) ,intent(in)    :: rhs
    allocate(lhs%Import_id,source=rhs)
  end subroutine

  integer(c_int) function NumSameIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumSameIDs=Epetra_Import_NumSameIDs(this%Import_id)
  end function

  integer(c_int) function NumPermuteIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumPermuteIDs=Epetra_Import_NumPermuteIDs(this%Import_id)
  end function

  integer(c_int) function NumRemoteIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumRemoteIDs=Epetra_Import_NumRemoteIDs(this%Import_id)
  end function

  integer(c_int) function NumExportIDs(this)
    class(Epetra_Import), intent(in) :: this
    NumExportIDs=Epetra_Import_NumExportIDs(this%Import_id)
  end function

  integer(c_int) function NumSend(this)
    class(Epetra_Import), intent(in) :: this
    NumSend=Epetra_Import_NumSend(this%Import_id)
  end function

  integer(c_int) function NumRecv(this)
    class(Epetra_Import), intent(in) :: this
    NumRecv=Epetra_Import_NumRecv(this%Import_id)
  end function

  type(Epetra_BlockMap) function SourceMap(this)
   class(Epetra_Import), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: SourceMap_id
   SourceMap_id=Epetra_Import_SourceMap(this%Import_id)
   SourceMap=Epetra_BlockMap(SourceMap_id)
  end function

  type(Epetra_BlockMap) function TargetMap(this)
   class(Epetra_Import), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: TargetMap_id
   TargetMap_id=Epetra_Import_TargetMap(this%Import_id)
   TargetMap=Epetra_BlockMap(TargetMap_id)
  end function

  subroutine finalize(this)
    type(Epetra_Import) :: this
    print *,'finalize_Import'
    call Epetra_Import_Destroy( this%Import_id ) 
    deallocate (this%Import_id)
  end subroutine

  subroutine force_finalization(this)
    class(Epetra_Import) ,intent(inout) :: this
    if (associated(this%Import_id)) then
      call finalize(this) 
    else
      print *,' finalization for Epetra_Import received  with unassociated object'
    end if
  end subroutine

end module 

