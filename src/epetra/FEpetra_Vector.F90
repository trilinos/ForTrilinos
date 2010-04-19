module FEpetra_Vector
  use ForTrilinos_enums   ,only: FT_Epetra_MultiVector_ID_t,FT_Epetra_Vector_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t,FT_boolean_t 
  use ForTrilinos_table_man
  use FEpetra_MultiVector ,only: Epetra_MultiVector
  use FEpetra_BlockMap    !,only: Epetra_BlockMap !use to circumvent reported compiler bug
  use iso_c_binding       ,only: c_int
  use forepetra
  private                                    ! Hide everything by default
  public :: Epetra_Vector,Epetra_MultiVector ! Expose type/constructors/methods
  implicit none

  type ,extends(Epetra_MultiVector)      :: Epetra_Vector !"shell"
    private
    type(FT_Epetra_Vector_ID_t) ,pointer :: vector_id => null()
  contains
     procedure         :: get_EpetraVector_ID 
     procedure ,nopass :: alias_EpetraVector_ID
     procedure         :: generalize 
     procedure         :: assign_to_Epetra_Vector
     generic :: assignment(=) => assign_to_Epetra_Vector
     procedure         :: force_finalization 
     final :: finalize
  end type

   interface Epetra_Vector ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface
 
contains

  type(FT_Epetra_Vector_ID_t) function from_struct(id)
     type(FT_Epetra_Vector_ID_t) ,intent(in) :: id
     from_struct = id
  end function

  ! Original C++ prototype:
  ! Epetra_Vector(const Epetra_BlockMap& Map, bool zeroOut = true);
  ! CTrilinos prototype:
  ! CT_Epetra_Vector_ID_t Epetra_Vector_Create ( CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut );

  type(FT_Epetra_Vector_ID_t) function from_scratch(BlockMap,zero_initial)
    use ForTrilinos_enums ,only: FT_boolean_t,FT_FALSE,FT_TRUE,FT_Epetra_BlockMap_ID_t
    use FEpetra_BlockMap  ,only: Epetra_BlockMap
    class(Epetra_BlockMap) ,intent(in) :: BlockMap
    logical ,optional      ,intent(in) :: zero_initial
    integer(FT_boolean_t)              :: zero_out
    if (present(zero_initial).and.zero_initial) then
     zero_out=FT_TRUE
    else
     zero_out=FT_FALSE
    endif
    from_scratch = Epetra_Vector_Create(BlockMap%get_EpetraBlockMap_ID(),zero_out)
  end function

  ! Original C++ prototype:
  ! Epetra_Vector(const Epetra_Vector& Source);
  ! CTrilinos prototype:
  ! CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( CT_Epetra_Vector_ID_t SourceID );

  type(FT_Epetra_Vector_ID_t) function duplicate(original)
    type(Epetra_Vector) ,intent(in) :: original
    duplicate = Epetra_Vector_Duplicate(original%vector_id)
  end function

  type(FT_Epetra_Vector_ID_t) function get_EpetraVector_ID(this)
    class(Epetra_Vector) ,intent(in) :: this 
    if (associated(this%vector_id)) then
     get_EpetraVector_ID=this%vector_id
    else
     stop 'get_EpetraVector_ID: vector_id is unassociated'
    end if
  end function
 
  type(FT_Epetra_Vector_ID_t) function alias_EpetraVector_ID(generic_id)
    use iso_c_binding        ,only: c_loc
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t, FT_Epetra_Vector_ID
    use ForTrilinos_table_man,only: CT_Alias 
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_Vector_ID))
    alias_EpetraVector_ID=degeneralize_EpetraVector(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(Epetra_Vector) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%vector_id))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_Vector) ,intent(in) ,target :: this
   ! generalize = Epetra_Vector_Generalize ( this%vector_id )
   ! ____ Use for CTrilinos function implementation ______
  end function
  
  type(FT_Epetra_Vector_ID_t) function degeneralize_EpetraVector(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums           ,only: ForTrilinos_Universal_ID_t,FT_Epetra_Vector_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                 ,value   :: generic_id
    type(FT_Epetra_Vector_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraVector = local_ptr
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
   !degeneralize_EpetraVector = Epetra_Vector_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine assign_to_Epetra_Vector(lhs,rhs)
    class(Epetra_Vector)        ,intent(inout) :: lhs
    type(FT_Epetra_Vector_ID_t) ,intent(in)    :: rhs
    allocate(lhs%vector_id,source=rhs)
    lhs%Epetra_MultiVector=Epetra_MultiVector(lhs%alias_EpetraMultiVector_ID(lhs%generalize()))
  end subroutine

  subroutine finalize(this)
    type(Epetra_Vector) :: this
    call Epetra_Vector_Destroy( this%vector_id ) 
    deallocate(this%vector_id)
  end subroutine

  subroutine force_finalization(this)
    class(Epetra_Vector) ,intent(inout) :: this
    call this%Epetra_MultiVector%force_finalization()
    if (associated(this%vector_id)) then
      call finalize(this) 
    else
      print *,' finalization for Epetra_Vector received object with unassociated vector_id'
    end if
  end subroutine

end module 

