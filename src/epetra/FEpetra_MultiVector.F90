!*********************************************************************
! ForTrilinos: Object-Oriented Fortran 2003 interface to Trilinos
!                Copyright 2010 Sandia Corporation
!
! Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
! the U.S. Government retains certain rights in this software.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright
!    notice, this list of conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright
!    notice, this list of conditions and the following disclaimer in the
!    documentation and/or other materials provided with the distribution.
!
! 3. Neither the name of the Corporation nor the names of the
!    contributors may be used to endorse or promote products derived from
!    this software without specific prior written permission.
!
! THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
! PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
! CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
! EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
! PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
! PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
! LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
! NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
! SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! Questions? Contact Karla Morris  (knmorri@sandia.gov) or 
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************


module FEpetra_MultiVector
  use ForTrilinos_enums ,only: FT_Epetra_MultiVector_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_BlockMap  ,only: Epetra_BlockMap
  use iso_c_binding     ,only: c_int,c_double,c_char
  use forepetra
  implicit none
  private                      ! Hide everything by default
  public :: Epetra_MultiVector ! Expose type/constructors/methods

  type ,extends(universal)                    :: Epetra_MultiVector !"shell"
    private
    type(FT_Epetra_MultiVector_ID_t) :: MultiVector_id 
  contains
     procedure         :: remote_dealloc_EpetraMultiVector
     procedure         :: remote_dealloc
     procedure         :: get_EpetraMultiVector_ID 
     procedure ,nopass :: alias_EpetraMultiVector_ID
     procedure         :: generalize 
     ! Post-construction modification procedure 
     procedure         :: ReplaceGlobalValue_GlobalRow
     procedure         :: ReplaceGlobalValue_GlobalBlockRow
     generic :: ReplaceGlobalValue=>ReplaceGlobalValue_GlobalRow!,ReplaceGlobalValue_GlobalBlockRow
     !Mathematical Methods
     procedure         :: Dot
     procedure         :: Abs
     procedure         :: Reciprocal
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
     procedure         :: MinValue
     procedure         :: MaxValue
     procedure         :: MeanValue
     procedure         :: Multiply_Matrix
     procedure         :: Multiply_ByEl
     generic :: Multiply => Multiply_Matrix,Multiply_ByEl
     procedure         :: ReciprocalMultiply
     !Attribute Access Functions
     procedure         :: NumVectors
     procedure         :: MyLength
     procedure         :: GlobalLength
     procedure         :: Stride
     procedure         :: ConstantStride
  end type

   interface Epetra_MultiVector ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface

contains
  type(Epetra_MultiVector) function from_struct(id)
     type(FT_Epetra_MultiVector_ID_t) ,intent(in) :: id
     from_struct%MultiVector_id = id
     call from_struct%register_self
  end function

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut = true);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut ); 

  type(Epetra_MultiVector) function from_scratch(BlockMap,Num_Vectors,zero)
   use ForTrilinos_enums ,only: FT_boolean_t,FT_TRUE,FT_FALSE
   use iso_c_binding     ,only: c_int
   class(Epetra_BlockMap) ,intent(in) :: BlockMap
   integer(c_int)         ,intent(in) :: Num_Vectors
   logical  ,intent(in) :: zero 
   integer(FT_boolean_t) :: zeroOut 
   type(FT_Epetra_MultiVector_ID_t) :: from_scratch_id
   if (zero) zeroOut=FT_TRUE
   if (.not.zero) zeroOut=FT_FALSE
   from_scratch_id = Epetra_MultiVector_Create(BlockMap%get_EpetraBlockMap_ID(),Num_Vectors,zeroOut)
   call from_scratch%register_self
  end function

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_MultiVector& Source);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( CT_Epetra_MultiVector_ID_t SourceID );

  type(Epetra_MultiVector) function duplicate(this)
    type(Epetra_MultiVector) ,intent(in) :: this
    type(FT_Epetra_MultiVector_ID_t) :: duplicate_id
    duplicate_id = Epetra_MultiVector_Duplicate(this%MultiVector_id)
    call duplicate%register_self
  end function

  type(FT_Epetra_MultiVector_ID_t) function get_EpetraMultiVector_ID(this)
    class(Epetra_MultiVector) ,intent(in) :: this 
    get_EpetraMultiVector_ID=this%MultiVector_id
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
   error_out=Epetra_MultiVector_Dot(this%MultiVector_id,x%MultiVector_id,result)
   if (present(err)) err=error(error_out)
  end subroutine

  subroutine Abs(this,A,err)
   class(Epetra_MultiVector), intent(in) :: this
   class(Epetra_MultiVector), intent(in) :: A 
   type(error),optional,intent(out)   :: err
   integer(c_int)                        :: error_out
   error_out=Epetra_MultiVector_Abs(this%MultiVector_id,A%MultiVector_id)
   if (present(err)) err=error(error_out)
  end subroutine

  subroutine Reciprocal(this,A,err)
   class(Epetra_MultiVector), intent(in) :: this
   class(Epetra_MultiVector), intent(in) :: A
   type(error),optional,intent(out)   :: err
   integer(c_int)                        :: error_out
   error_out=Epetra_MultiVector_Reciprocal(this%MultiVector_id,A%MultiVector_id)
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

  function MinValue(this,err) result(MinValue_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: MinValue_val
    integer(c_int)                           :: error_out
    allocate(MinValue_val(this%NumVectors()))
    error_out = Epetra_MultiVector_MinValue(this%MultiVector_id,MinValue_val)
    if (present(err)) err=error(error_out)
  end function
 
 function MaxValue(this,err) result(MaxValue_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: MaxValue_val
    integer(c_int)                           :: error_out
    allocate(MaxValue_val(this%NumVectors()))
    error_out = Epetra_MultiVector_MaxValue(this%MultiVector_id,MaxValue_val)
    if (present(err)) err=error(error_out)
  end function

 function MeanValue(this,err) result(MeanValue_val)
    class(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: MeanValue_val
    integer(c_int)                           :: error_out
    allocate(MeanValue_val(this%NumVectors()))
    error_out = Epetra_MultiVector_MeanValue(this%MultiVector_id,MeanValue_val)
    if (present(err)) err=error(error_out)
  end function

  subroutine Multiply_Matrix(this,TransA,TransB,ScalarAB,A,B,ScalarThis,err) 
    class(Epetra_MultiVector)   ,intent(in) :: this
    character(c_char), intent(in) :: TransA
    character(c_char), intent(in) :: TransB
    class(Epetra_MultiVector)   ,intent(in)  :: A,B
    real(c_double), intent(in) :: ScalarAB, ScalarThis 
    type(error) ,optional    ,intent(out) :: err
    integer(c_int)                           :: error_out
    error_out = Epetra_MultiVector_Multiply_Matrix(this%MultiVector_id,TransA,TransB,ScalarAB,&
                                 A%MultiVector_id,B%MultiVector_id,ScalarThis)
    if (present(err)) err=error(error_out)
  end subroutine

 subroutine Multiply_ByEl(this,ScalarAB,A,B,ScalarThis,err)
    class(Epetra_MultiVector)   ,intent(in) :: this
    class(Epetra_MultiVector)   ,intent(in)  :: A,B
    real(c_double), intent(in) :: ScalarAB, ScalarThis
    type(error) ,optional    ,intent(out) :: err
    integer(c_int)                           :: error_out
    error_out = Epetra_MultiVector_Multiply_ByEl(this%MultiVector_id,ScalarAB,&
                                 A%MultiVector_id,B%MultiVector_id,ScalarThis)
    if (present(err)) err=error(error_out)
  end subroutine

  subroutine ReciprocalMultiply(this,ScalarAB,A,B,ScalarThis,err)
    class(Epetra_MultiVector)   ,intent(in) :: this
    class(Epetra_MultiVector)   ,intent(in)  :: A,B
    real(c_double), intent(in) :: ScalarAB, ScalarThis
    type(error) ,optional    ,intent(out) :: err
    integer(c_int)                           :: error_out
    error_out = Epetra_MultiVector_ReciprocalMultiply(this%MultiVector_id,ScalarAB,&
                                 A%MultiVector_id,B%MultiVector_id,ScalarThis)
    if (present(err)) err=error(error_out)
  end subroutine

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

  subroutine remote_dealloc_EpetraMultiVector(this)
    class(Epetra_MultiVector),intent(inout) :: this
    call Epetra_MultiVector_Destroy( this%MultiVector_id ) 
  end subroutine

  subroutine remote_dealloc(this)
    class(Epetra_MultiVector),intent(inout) :: this
    call Epetra_MultiVector_Destroy( this%MultiVector_id ) 
  end subroutine

end module 

