module FEpetra_Export
  use ForTrilinos_enums ,only: FT_Epetra_Comm_ID_t,FT_Epetra_Export_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal
  use FEpetra_Comm  ,only: Epetra_Comm
  use FEpetra_BlockMap ,only: Epetra_BlockMap
  use iso_c_binding ,only: c_int
  use forepetra
  implicit none
  private                   ! Hide everything by default
  public :: Epetra_Export ! Expose type/constructors/methods

  type ,extends(universal)                 :: Epetra_Export !"shell"
    private
    type(FT_Epetra_Export_ID_t) ,pointer :: Export_id => null()
  contains
     !Developers only
     procedure         :: get_EpetraExport_ID 
     procedure ,nopass :: alias_EpetraExport_ID
     procedure         :: generalize 
     procedure         :: assign_to_Epetra_Export
     generic :: assignment(=) => assign_to_Epetra_Export
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

   interface Epetra_Export ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface
 
contains
  type(FT_Epetra_Export_ID_t) function from_struct(id)
     type(FT_Epetra_Export_ID_t) ,intent(in) :: id
     from_struct = id
  end function
 
  ! Original C++ prototype:
  ! Epetra_Export( const Epetra_BlockMap & SourceMap, const Epetra_BlockMap & TargetMap );
  ! CTrilinos prototype:
  ! CT_Epetra_Export_ID_t Epetra_Export_Create ( CT_Epetra_BlockMap_ID_t SourceMapID, CT_Epetra_BlockMap_ID_t TargetMapID );

  type(FT_Epetra_Export_ID_t) function from_scratch(SourceMap,TargetMap)
   !use ForTrilinos_enums ,only : FT_Epetra_Comm_ID_t,FT_Epetra_BlockMap_ID_t
    type(Epetra_BlockMap), intent(in) :: TargetMap
    type(Epetra_BlockMap), intent(in) :: SourceMap
    from_scratch = Epetra_Export_Create(SourceMap%get_EpetraBlockMap_ID(),TargetMap%get_EpetraBlockMap_ID())
  end function

  ! Original C++ prototype:
  ! Epetra_Export(const Epetra_Import& Importer);
  ! CTrilinos prototype:
  ! CT_Epetra_Export_ID_t Epetra_Export_Duplicate ( CT_Epetra_Import_ID_t ImporterID );

  type(FT_Epetra_Export_ID_t) function duplicate(original)
    type(Epetra_Export) ,intent(in) :: original
    duplicate = Epetra_Export_Duplicate(original%Export_id)
  end function

  type(FT_Epetra_Export_ID_t) function get_EpetraExport_ID(this)
    class(Epetra_Export) ,intent(in) :: this 
    if (associated(this%Export_id)) then
     get_EpetraExport_ID=this%Export_id
    else
     stop 'get_EpetraExport_ID: Export_id is unassociated'
    end if
  end function
  
  type(FT_Epetra_Export_ID_t) function alias_EpetraExport_ID(generic_id)
    use ForTrilinos_table_man
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Epetra_Export_ID
    use iso_c_binding     ,only: c_loc
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Export_ID))
    alias_EpetraExport_ID=degeneralize_EpetraExport(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(Epetra_Export) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%Export_ID))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(Epetra_Export) ,intent(in) ,target :: this
   !generalize = Epetra_Export_Generalize ( this%Export_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_Export_ID_t) function degeneralize_EpetraExport(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_Export_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                   ,value   :: generic_id
    type(FT_Epetra_Export_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraExport = local_ptr
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
   !degeneralize_EpetraExport = Epetra_Export_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  subroutine assign_to_Epetra_Export(lhs,rhs)
    class(Epetra_Export)        ,intent(inout) :: lhs
    type(FT_Epetra_Export_ID_t) ,intent(in)    :: rhs
    allocate(lhs%Export_id,source=rhs)
  end subroutine

  integer(c_int) function NumSameIDs(this)
    class(Epetra_Export), intent(in) :: this
    NumSameIDs=Epetra_Export_NumSameIDs(this%Export_id)
  end function

  integer(c_int) function NumPermuteIDs(this)
    class(Epetra_Export), intent(in) :: this
    NumPermuteIDs=Epetra_Export_NumPermuteIDs(this%Export_id)
  end function

  integer(c_int) function NumRemoteIDs(this)
    class(Epetra_Export), intent(in) :: this
    NumRemoteIDs=Epetra_Export_NumRemoteIDs(this%Export_id)
  end function

  integer(c_int) function NumExportIDs(this)
    class(Epetra_Export), intent(in) :: this
    NumExportIDs=Epetra_Export_NumExportIDs(this%Export_id)
  end function

  integer(c_int) function NumSend(this)
    class(Epetra_Export), intent(in) :: this
    NumSend=Epetra_Export_NumSend(this%Export_id)
  end function

  integer(c_int) function NumRecv(this)
    class(Epetra_Export), intent(in) :: this
    NumRecv=Epetra_Export_NumRecv(this%Export_id)
  end function

  type(Epetra_BlockMap) function SourceMap(this)
   class(Epetra_Export), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: SourceMap_id
   SourceMap_id=Epetra_Export_SourceMap(this%Export_id)
   SourceMap=Epetra_BlockMap(SourceMap_id)
  end function

  type(Epetra_BlockMap) function TargetMap(this)
   class(Epetra_Export), intent(in) :: this
   type(FT_Epetra_BlockMap_ID_t) :: TargetMap_id
   TargetMap_id=Epetra_Export_TargetMap(this%Export_id)
   TargetMap=Epetra_BlockMap(TargetMap_id)
  end function

  subroutine finalize(this)
    type(Epetra_Export) :: this
    print *,'finalize_Export'
    call Epetra_Export_Destroy( this%Export_id ) 
    deallocate (this%Export_id)
  end subroutine

  subroutine force_finalization(this)
    class(Epetra_Export) ,intent(inout) :: this
    if (associated(this%Export_id)) then
      call finalize(this) 
    else
      print *,' finalization for Epetra_Export received  with unassociated object'
    end if
  end subroutine

end module 

