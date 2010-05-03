module FEpetra_MultiVector
  use ForTrilinos_enums ,only: FT_Epetra_MultiVector_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal
  use ForTrilinos_error
  use FEpetra_BlockMap  ,only: Epetra_BlockMap
  use iso_c_binding     ,only: c_int,c_double,c_char
  use forepetra
  implicit none
  private                      ! Hide everything by default
  public :: Epetra_MultiVector ! Expose type/constructors/methods

  type ,extends(universal)                    :: Epetra_MultiVector !"shell"
    private
    type(FT_Epetra_MultiVector_ID_t) ,pointer :: MultiVector_id => null()
  contains
     procedure         :: get_EpetraMultiVector_ID 
     procedure ,nopass :: alias_EpetraMultiVector_ID
     procedure         :: generalize 
     procedure         :: assign_to_Epetra_MultiVector
     generic :: assignment(=) => assign_to_Epetra_MultiVector
     ! Post-construction modification procedure 
     procedure         :: ReplaceGlobalValue_GlobalRow
     procedure         :: ReplaceGlobalValue_GlobalBlockRow
     generic :: ReplaceGlobalValue=>ReplaceGlobalValue_GlobalRow!,ReplaceGlobalValue_GlobalBlockRow
     !Mathematical Methods
     procedure         :: Dot
     procedure         :: Scale_Self
     procedure         :: Scale_Other
     generic :: Scale => Scale_Self,Scale_Other
     procedure         :: PutScalar
     procedure         :: Random
     ! Extraction Methods
     procedure         :: ExtractCopy_2DA
     generic :: ExtractCopy=>ExtractCopy_2DA
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

   interface Epetra_MultiVector ! constructors
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

  type(FT_Epetra_MultiVector_ID_t) function from_scratch(BlockMap,Num_Vectors,zero)
   use ForTrilinos_enums ,only: FT_boolean_t,FT_TRUE,FT_FALSE
   use iso_c_binding     ,only: c_int
   class(Epetra_BlockMap) ,intent(in) :: BlockMap
   integer(c_int)         ,intent(in) :: Num_Vectors
   logical  ,intent(in) :: zero 
   integer(FT_boolean_t) :: zeroOut 
   if (zero) zeroOut=FT_TRUE
   if (.not.zero) zeroOut=FT_FALSE
   from_scratch = Epetra_MultiVector_Create(BlockMap%get_EpetraBlockMap_ID(),Num_Vectors,zeroOut)
  end function

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_MultiVector& Source);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( CT_Epetra_MultiVector_ID_t SourceID );

  type(FT_Epetra_MultiVector_ID_t) function duplicate(original)
    type(Epetra_MultiVector) ,intent(in) :: original
    duplicate = Epetra_MultiVector_Duplicate(original%MultiVector_id)
  end function

  type(FT_Epetra_MultiVector_ID_t) function get_EpetraMultiVector_ID(this)
    class(Epetra_MultiVector) ,intent(in) :: this 
    if (associated(this%MultiVector_id)) then
     get_EpetraMultiVector_ID=this%MultiVector_id
    else
     stop 'get_EpetraMultiVector_ID: MultiVector_id is unassociated'
    end if
  end function
  
  type(FT_Epetra_MultiVector_ID_t) function alias_EpetraMultiVector_ID(generic_id)
    use ForTrilinos_table_man,only: CT_Alias
    use iso_c_binding        ,only: c_loc
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t,FT_Epetra_MultiVector_ID
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
   class(Epetra_MultiVector) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%MultiVector_ID))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_MultiVector) ,intent(in) ,target :: this
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
 
  subroutine assign_to_Epetra_MultiVector(lhs,rhs)
    class(Epetra_MultiVector)       ,intent(inout):: lhs
    type(FT_Epetra_MultiVector_ID_t),intent(in)   :: rhs
    allocate(lhs%MultiVector_id,source=rhs)
  end subroutine

  subroutine ReplaceGlobalValue_GlobalRow(this,GlobalRow,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalRow
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_mod ! To account for Fortran index base 1
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_mod = VectorIndex-1
    error_out=Epetra_MultiVector_ReplaceGlobalValue(this%MultiVector_id,GlobalRow,VectorIndex_mod,ScalarValue)
    if (present(err)) err=error(error_out)
  end subroutine
  
  subroutine ReplaceGlobalValue_GlobalBlockRow(this,GlobalRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    class(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalRow
    integer(c_int), intent(in) :: BlockRowOffset
    integer(c_int), intent(in) :: VectorIndex 
    integer(c_int)             :: VectorIndex_mod ! To account for Fortran index base 1
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
    integer(c_int)     :: error_out
    VectorIndex_mod = VectorIndex-1
    error_out=Epetra_MultiVector_ReplaceGlobalValue_BlockPos(this%MultiVector_id,GlobalRow,BlockRowOffset,VectorIndex_mod,ScalarValue)
    if (present(err)) err=error(error_out)
  end subroutine

  subroutine Dot(this,x,result,err)
   class(Epetra_MultiVector), intent(in) :: this
   class(Epetra_MultiVector), intent(in) :: x
   real(c_double),dimension(:)           :: result
   type(error),optional,intent(out)   :: err
   integer(c_int)                        :: error_out
   error_out=Epetra_multiVector_Dot(this%MultiVector_id,x%MultiVector_id,result)
   if (present(err)) err=error(error_out)
  end subroutine
 
  subroutine Scale_Self(this,scalar_value,err)
    class(Epetra_MultiVector), intent(inout) :: this
    real(c_double),            intent(in)    :: scalar_value
    type(error),optional,intent(out)      :: err
    integer(c_int)                           :: error_out
    error_out=Epetra_MultiVector_Scale_Self(this%MultiVector_id,scalar_value)
    if (present(err)) err=error(error_out)
  end subroutine
  
  subroutine Scale_Other(this,scalar_value,MultiVector,err)
    class(Epetra_MultiVector), intent(inout) :: this
    class(Epetra_MultiVector), intent(in)    :: MultiVector 
    real(c_double),            intent(in)    :: scalar_value
    type(error),optional,intent(out)      :: err
    integer(c_int)                           :: error_out
    error_out=Epetra_MultiVector_Scale(this%MultiVector_id,scalar_value,MultiVector%MultiVector_id)
    if (present(err)) err=error(error_out)
  end subroutine

  subroutine PutScalar(this,scalar,err)
   class(Epetra_MultiVector) ,intent(inout) :: this
   type(error) ,optional  ,intent(out)   :: err
   real(c_double)            ,intent(in)    :: scalar
   integer(c_int)                           :: error_out
   error_out=Epetra_MultiVector_PutScalar(this%MultiVector_id,scalar)
   if (present(err)) err=error(error_out)
  end subroutine
 
  subroutine Random(this,err)
   class(Epetra_MultiVector) ,intent(inout) :: this
   type(error) ,optional  ,intent(out)   :: err
   integer(c_int)                           :: error_out
   error_out = Epetra_MultiVector_Random (this%MultiVector_id)
   if (present(err)) err=error(error_out)
  end subroutine

  function ExtractCopy_2DA(this,MyLDA,err) result(A)
   class(Epetra_MultiVector),intent(in) :: this
   real(c_double),dimension(:,:),allocatable :: A
   integer(c_int),intent(in) :: MyLDA
   type(error),optional,intent(out) :: err
   integer(c_int)           :: error_out
   allocate(A(this%NumVectors(),this%MyLength()))
   error_out=Epetra_MultiVector_ExtractCopy_Fill2DA(this%MultiVector_id,A,MyLDA)
   if (present(err)) err=error(error_out)
  end function

  subroutine Update_WithA(this,scalarA,A,scalarThis,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    real(c_double)            ,intent(in) :: scalarA
    class(Epetra_MultiVector) ,intent(in) :: A
    real(c_double)            ,intent(in) :: scalarThis
    type(error) ,optional  ,intent(out):: err
    integer(c_int)                        :: error_out
    error_out = Epetra_MultiVector_Update_WithA(this%MultiVector_id,scalarA,A%MultiVector_id,scalarThis)
     if (present(err)) err=error(error_out)
  end subroutine 

  subroutine Update_WithAB(this,scalarA,A,scalarB,B,scalarThis,err)
    class(Epetra_MultiVector) ,intent(inout) :: this
    real(c_double)            ,intent(in) :: scalarA
    class(Epetra_MultiVector) ,intent(in) :: A
    real(c_double)            ,intent(in) :: scalarB
    class(Epetra_MultiVector) ,intent(in) :: B
    real(c_double)            ,intent(in) :: scalarThis
    type(error) ,optional  ,intent(out):: err
    integer(c_int)                        :: error_out
    error_out = Epetra_MultiVector_Update_WithAB(this%MultiVector_id,scalarA,A%MultiVector_id,scalarB,B%MultiVector_id,scalarThis)
    if (present(err)) err=error(error_out)
  end subroutine

  function Norm1(this,err) result(Norm1_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: Norm1_val 
    integer(c_int)                           :: error_out
    allocate(Norm1_val(this%NumVectors()))
    error_out = Epetra_MultiVector_Norm1(this%MultiVector_id,Norm1_val)
    if (present(err)) err=error(error_out)
  end function 
 
  function Norm2(this,err) result(Norm2_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: Norm2_val 
    integer(c_int)                           :: error_out
    allocate(Norm2_val(this%NumVectors()))
    error_out = Epetra_MultiVector_Norm2(this%MultiVector_id,Norm2_val)
    if (present(err)) err=error(error_out)
  end function
 
  function NormInf(this,err) result(NormInf_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: NormInf_val 
    integer(c_int)                           :: error_out
    allocate(NormInf_val(this%NumVectors()))
    error_out = Epetra_MultiVector_NormInf(this%MultiVector_id,NormInf_val)
    if (present(err)) err=error(error_out)
  end function 

  function NormWeighted(this,weights,err) result(NormWeighted_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    class(Epetra_MultiVector)   ,intent(in)  :: weights 
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: NormWeighted_val 
    integer(c_int)                           :: error_out
    allocate(NormWeighted_val(this%NumVectors()))
    error_out = Epetra_MultiVector_NormWeighted(this%MultiVector_id,weights%MultiVector_id,NormWeighted_val)
    if (present(err)) err=error(error_out)
  end function 

  integer(c_int) function NumVectors(this)
    class(Epetra_MultiVector) ,intent(in) :: this
    NumVectors=Epetra_MultiVector_NumVectors(this%MultiVector_id)
  end function 

  integer(c_int) function MyLength(this)
    class(Epetra_MultiVector) ,intent(in) :: this
    MyLength=Epetra_MultiVector_MyLength(this%MultiVector_id)
  end function 

  integer(c_int) function GlobalLength(this)
    class(Epetra_MultiVector) ,intent(in) :: this
    GlobalLength=Epetra_MultiVector_GlobalLength(this%MultiVector_id)
  end function 

  integer(c_int) function Stride(this)
    class(Epetra_MultiVector) ,intent(in) :: this
    Stride=Epetra_MultiVector_Stride(this%MultiVector_id)
  end function 

  integer(FT_boolean_t) function ConstantStride(this)
    use ForTrilinos_enums ,only:FT_Epetra_MultiVector_ID_t,FT_boolean_t
    class(Epetra_MultiVector) ,intent(in) :: this
    ConstantStride=Epetra_MultiVector_ConstantStride(this%MultiVector_id)
  end function 

  subroutine finalize(this)
    type(Epetra_MultiVector) :: this
    call Epetra_MultiVector_Destroy( this%MultiVector_id ) 
    deallocate (this%MultiVector_id)
  end subroutine

  subroutine force_finalization(this)
    class(Epetra_MultiVector) ,intent(inout) :: this
    if (associated(this%MultiVector_id)) then
      call finalize(this) 
    else
      print *,' finalization for Epetra_MultiVector received  with unassociated object'
    end if
  end subroutine

end module 

