module FEpetra_MultiVector
  use ForTrilinos_enums ,only: FT_Epetra_MultiVector_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_universal
  use FEpetra_BlockMap  ,only: epetra_BlockMap
  use iso_c_binding     ,only: c_int,c_double
  use forepetra
  implicit none
  private                      ! Hide everything by default
  public :: epetra_MultiVector ! Expose type/constructors/methods

  type ,extends(universal)                    :: epetra_MultiVector !"shell"
    private
    type(FT_Epetra_MultiVector_ID_t) ,pointer :: MultiVector_id => null()
  contains
     procedure         :: get_EpetraMultiVector_ID 
     procedure ,nopass :: alias_EpetraMultiVector_ID
     procedure         :: generalize 
     procedure         :: assign_to_epetra_MultiVector
     generic :: assignment(=) => assign_to_epetra_MultiVector
     ! Post-construction modification procedure 
     procedure         :: Random
     ! Extraction Methods
     ! Mathematical Methods
     procedure         :: Update_WithA
     procedure         :: Update_WithAB
     generic :: Update => Update_WithA,Update_WithAB
     procedure         :: Norm1
     procedure         :: Norm2
     procedure         :: NormInf
     procedure         :: NormWeighted
     !Attribute Access Functions
     procedure         :: NumVectors
     procedure         :: MyLength
     procedure         :: GlobalLength
     procedure         :: Stride
     procedure         :: ConstantStride
     !Memory Management
     procedure         :: force_finalization 
     final :: finalize
  end type

   interface epetra_MultiVector ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface

contains
  type(FT_Epetra_MultiVector_ID_t) function from_struct(id)
     type(FT_Epetra_MultiVector_ID_t) ,intent(in) :: id
     from_struct = id
  end function

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut = true);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut ); 

  type(FT_Epetra_MultiVector_ID_t) function from_scratch(BlockMap,Num_Vectors,zeroOut)
   use ForTrilinos_enums ,only: FT_boolean_t
   use iso_c_binding     ,only: c_int
   class(epetra_BlockMap) ,intent(in) :: BlockMap
   integer(c_int)         ,intent(in) :: Num_Vectors
   integer(FT_boolean_t)  ,intent(in) :: zeroOut 
   from_scratch = Epetra_MultiVector_Create(BlockMap%get_EpetraBlockMap_ID(),Num_Vectors,zeroOut)
  end function

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_MultiVector& Source);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( CT_Epetra_MultiVector_ID_t SourceID );

  type(FT_Epetra_MultiVector_ID_t) function duplicate(original)
    type(epetra_MultiVector) ,intent(in) :: original
    duplicate = Epetra_MultiVector_Duplicate(original%MultiVector_id)
  end function

  type(FT_Epetra_MultiVector_ID_t) function get_EpetraMultiVector_ID(this)
    class(epetra_MultiVector) ,intent(in) :: this 
    if (associated(this%MultiVector_id)) then
     get_EpetraMultiVector_ID=this%MultiVector_id
    else
     stop 'get_EpetraMultiVector_ID: MultiVector_id is unassociated'
    end if
  end function
  
  type(FT_Epetra_MultiVector_ID_t) function alias_EpetraMultiVector_ID(generic_id)
    use ForTrilinos_utils ,only: CT_Alias
    use iso_c_binding     ,only: c_loc
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Epetra_MultiVector_ID
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_MultiVector_ID))
    alias_EpetraMultiVector_ID=degeneralize_EpetraMultiVector(c_loc(alias_id))
    deallocate(alias_id)
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only : c_loc
   class(epetra_MultiVector) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%MultiVector_ID))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(epetra_MultiVector) ,intent(in) ,target :: this
   ! generalize = Epetra_MultiVector_Generalize ( this%MultiVector_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Epetra_MultiVector_ID_t) function degeneralize_EpetraMultiVector(generic_id) bind(C)
  ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_MultiVector_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                      ,value   :: generic_id
    type(FT_Epetra_MultiVector_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraMultiVector = local_ptr
  ! ____ Use for ForTrilinos function implementation ______

  ! ____ Use for CTrilinos function implementation ______
  !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
  !degeneralize_EpetraMultiVector = Epetra_MultiVector_Degeneralize(generic_id)
  ! ____ Use for CTrilinos function implementation ______
  end function
 
  subroutine assign_to_epetra_MultiVector(lhs,rhs)
    class(epetra_MultiVector), intent(out) :: lhs
    type(FT_Epetra_MultiVector_ID_t), intent(in) :: rhs
    allocate(lhs%MultiVector_id,source=rhs)
  end subroutine
  
  subroutine Random(this,error)
   class(epetra_MultiVector) ,intent(inout) :: this
   integer(c_int),intent(out),optional :: error
   integer(c_int)                      :: error_out
   error_out = Epetra_MultiVector_Random (this%MultiVector_id)
   if(present(error)) error=error_out
  end subroutine

  subroutine Update_WithA(this,scalarA,A,scalarThis,error)
    class(epetra_MultiVector) ,intent(inout) :: this
    real(c_double)            ,intent(in) :: scalarA
    class(epetra_MultiVector) ,intent(in) :: A
    real(c_double)            ,intent(in) :: scalarThis
    integer(c_int) ,optional  ,intent(out):: error
    integer(c_int)                        :: error_out
    error_out = Epetra_MultiVector_Update_WithA(this%MultiVector_id,scalarA,A%MultiVector_id,scalarThis)   
    if (present(error)) error=error_out
  end subroutine 

  subroutine Update_WithAB(this,scalarA,A,scalarB,B,scalarThis,error)
    class(epetra_MultiVector) ,intent(inout) :: this
    real(c_double)            ,intent(in) :: scalarA
    class(epetra_MultiVector) ,intent(in) :: A
    real(c_double)            ,intent(in) :: scalarB
    class(epetra_MultiVector) ,intent(in) :: B
    real(c_double)            ,intent(in) :: scalarThis
    integer(c_int) ,optional  ,intent(out):: error
    integer(c_int)                        :: error_out
    error_out = Epetra_MultiVector_Update_WithAB(this%MultiVector_id,scalarA,A%MultiVector_id,scalarB,B%MultiVector_id,scalarThis)
    if (present(error)) error=error_out
  end subroutine

  function Norm1(this,error) result(Norm1_val)
    class(epetra_MultiVector)   ,intent(in)  :: this
    integer(c_int) ,optional    ,intent(out) :: error
    real(c_double) ,dimension(:),allocatable :: Norm1_val 
    integer(c_int)                           :: error_out
    allocate(Norm1_val(this%NumVectors()))
    error_out = Epetra_MultiVector_Norm1(this%MultiVector_id,Norm1_val)
    if (present(error)) error=error_out
  end function 
 
  function Norm2(this,error) result(Norm2_val)
    class(epetra_MultiVector)   ,intent(in)  :: this
    integer(c_int) ,optional    ,intent(out) :: error
    real(c_double) ,dimension(:),allocatable :: Norm2_val 
    integer(c_int)                           :: error_out
    allocate(Norm2_val(this%NumVectors()))
    error_out = Epetra_MultiVector_Norm2(this%MultiVector_id,Norm2_val)
    if (present(error)) error=error_out
  end function
 
  function NormInf(this,error) result(NormInf_val)
    class(epetra_MultiVector)   ,intent(in)  :: this
    integer(c_int) ,optional    ,intent(out) :: error
    real(c_double) ,dimension(:),allocatable :: NormInf_val 
    integer(c_int)                           :: error_out
    allocate(NormInf_val(this%NumVectors()))
    error_out = Epetra_MultiVector_NormInf(this%MultiVector_id,NormInf_val)
    if (present(error)) error=error_out
  end function 

  function NormWeighted(this,weights,error) result(NormWeighted_val)
    class(epetra_MultiVector)   ,intent(in)  :: this
    class(epetra_MultiVector)   ,intent(in)  :: weights 
    integer(c_int) ,optional    ,intent(out) :: error
    real(c_double) ,dimension(:),allocatable :: NormWeighted_val 
    integer(c_int)                           :: error_out
    allocate(NormWeighted_val(this%NumVectors()))
    error_out = Epetra_MultiVector_NormWeighted(this%MultiVector_id,weights%MultiVector_id,NormWeighted_val)
    if (present(error)) error=error_out
  end function 

  integer(c_int) function NumVectors(this)
    class(epetra_MultiVector) ,intent(in) :: this
    NumVectors=Epetra_MultiVector_NumVectors(this%MultiVector_id)
  end function 

  integer(c_int) function MyLength(this)
    class(epetra_MultiVector) ,intent(in) :: this
    MyLength=Epetra_MultiVector_MyLength(this%MultiVector_id)
  end function 

  integer(c_int) function GlobalLength(this)
    class(epetra_MultiVector) ,intent(in) :: this
    GlobalLength=Epetra_MultiVector_GlobalLength(this%MultiVector_id)
  end function 

  integer(c_int) function Stride(this)
    class(epetra_MultiVector) ,intent(in) :: this
    Stride=Epetra_MultiVector_Stride(this%MultiVector_id)
  end function 

  integer(FT_boolean_t) function ConstantStride(this)
    use ForTrilinos_enums ,only:FT_Epetra_MultiVector_ID_t,FT_boolean_t
    class(epetra_MultiVector) ,intent(in) :: this
    ConstantStride=Epetra_MultiVector_ConstantStride(this%MultiVector_id)
  end function 

  subroutine finalize(this)
    type(epetra_MultiVector) :: this
    call Epetra_MultiVector_Destroy( this%MultiVector_id ) 
    deallocate (this%MultiVector_id)
  end subroutine

  subroutine force_finalization(this)
    class(epetra_MultiVector) ,intent(inout) :: this
    if (associated(this%MultiVector_id)) then
      call finalize(this) 
      deallocate(this%MultiVector_id)
    else
      print *,' finalization for epetra_MultiVector received  with unassociated object'
    end if
  end subroutine

end module 

