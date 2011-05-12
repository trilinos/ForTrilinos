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
  use ForTrilinos_enums ! ,only: FT_Epetra_MultiVector_ID_t,FT_Epetra_Map_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal,only:universal
  use ForTrilinos_error
  use FEpetra_DistObject, only: Epetra_DistObject
  use FEpetra_BlockMap  ,only: Epetra_BlockMap
  use iso_c_binding     ,only: c_int,c_double,c_char
  use forepetra
  implicit none
  private                      ! Hide everything by default
  public :: Epetra_MultiVector ! Expose type/constructors/methods

  type  Epetra_MultiVector !,extends(universal) "shell"
  contains
  end type

contains

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_BlockMap& Map, int NumVectors, bool zeroOut = true);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Create ( CT_Epetra_BlockMap_ID_t MapID, int NumVectors, boolean zeroOut ); 

  !> @name Constructor Functions
  !! @{
 
  !> <BR> Epetra_MultiVector constructor conformal to a BlockMap, optionally zero the newly created vector.
  type(Epetra_MultiVector) function Epetra_MultiVector(BlockMap,Num_Vectors,zero)
    type(Epetra_BlockMap) ,intent(in) :: BlockMap
    integer(c_int)         ,intent(in) :: Num_Vectors
    logical  ,intent(in) :: zero 
  end function

  ! Original C++ prototype:
  ! Epetra_MultiVector(const Epetra_MultiVector& Source);
  ! CTrilinos prototype:
  ! CT_Epetra_MultiVector_ID_t Epetra_MultiVector_Duplicate ( CT_Epetra_MultiVector_ID_t SourceID );

  !> @name Constructor Functions
  !! @{
 
  !> <BR> Epetra_MultiVector copy constructor.
  type(Epetra_MultiVector) function Epetra_MultiVector(this)
    type(Epetra_MultiVector) ,intent(in) :: this
  end function


  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces value at location (GlobalRow,VectorIndex) with ScalarValue
  subroutine ReplaceGlobalValue(this,GlobalRow,VectorIndex,ScalarValue,err)
    type(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalRow
    integer(c_int), intent(in) :: VectorIndex 
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Replaces value at location (GlobalBlockRow,BlockRowOffset,VectorIndex) with ScalarValue
  subroutine ReplaceGlobalValue(this,GlobalBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    type(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalBlockRow
    integer(c_int), intent(in) :: BlockRowOffset
    integer(c_int), intent(in) :: VectorIndex 
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
  end subroutine
 
  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Replaces value at location (MyRow,VectorIndex) with ScalarValue
  subroutine ReplaceMyValue(this,MyRow,VectorIndex,ScalarValue,err)
    type(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyRow
    integer(c_int), intent(in) :: VectorIndex 
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces value at location (MyBlockRow,BlockRowOffset,VectorIndex) with ScalarValue
  subroutine ReplaceMyValue(this,MyBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    type(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyBlockRow
    integer(c_int), intent(in) :: BlockRowOffset
    integer(c_int), intent(in) :: VectorIndex 
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
  end subroutine


  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Adds ScalarValue to value at location (GlobalRow,VectorIndex) 
  subroutine SumIntoGlobalValue ( this,GlobalRow,VectorIndex,ScalarValue,err)
    type(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalRow
    integer(c_int), intent(in) :: VectorIndex 
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Adds ScalarValue to value at location (GlobalBlockRow,BlockRowOffset,VectorIndex) 
  subroutine SumIntoGlobalValue ( this,GlobalBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    type(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: GlobalBlockRow
    integer(c_int), intent(in) :: BlockRowOffset
    integer(c_int), intent(in) :: VectorIndex 
    real(c_double), intent(in) :: ScalarValue 
  end subroutine

  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Adds ScalarValue to value at location (MyRow,VectorIndex) 
  subroutine SumIntoMyValue ( this,MyRow,VectorIndex,ScalarValue,err)
    type(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyRow
    integer(c_int), intent(in) :: VectorIndex 
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
  end subroutine
  
  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Adds ScalarValue to value at location (MyBlockRow,BlockRowOffset,VectorIndex) 
  subroutine SumIntoMyValue ( this,MyBlockRow,BlockRowOffset,VectorIndex,ScalarValue,err)
    type(Epetra_MultiVector), intent(in) :: this
    integer(c_int), intent(in) :: MyBlockRow
    integer(c_int), intent(in) :: BlockRowOffset
    integer(c_int), intent(in) :: VectorIndex 
    real(c_double), intent(in) :: ScalarValue 
    type(error),optional,intent(out) :: err
  end subroutine


  !> @name Post-construction modification routines
  !! @{

  !> <BR>  Replaces all entries with scalar
  subroutine PutScalar(this,scalar,err)
    type(Epetra_MultiVector) ,intent(inout) :: this
    type(error) ,optional  ,intent(out)   :: err
    real(c_double)            ,intent(in)    :: scalar
  end subroutine
 
  !> @name Post-construction modification routines
  !! @{

  !> <BR> Replaces all entries with random values
  subroutine Random(this,err)
    type(Epetra_MultiVector) ,intent(inout) :: this
    type(error) ,optional  ,intent(out)   :: err
  end subroutine



  !> @name Mathematical methods
  !! @{

  !> <BR> Computes the scalar product of corresponding pairs of vectors.
  function Dot(this,x,err) result(dot_)
    type(Epetra_MultiVector), intent(in) :: this
    type(Epetra_MultiVector), intent(in) :: x
    real(c_double),dimension(:),allocatable :: dot_ 
    type(error),optional,intent(out)   :: err
  end function 

  !> @name Mathematical methods
  !! @{

  !> <BR> Replaces target with element-wise absolute value of input
  subroutine Abs(this,A,err)
    type(Epetra_MultiVector), intent(inout) :: this
    type(Epetra_MultiVector), intent(in) :: A 
    type(error),optional,intent(out)   :: err
  end subroutine

  !> @name Mathematical methods
  !! @{

  !> <BR> Reciprocal  replaces target with element-wise reciprocal value of input
  subroutine Reciprocal(this,A,err)
    type(Epetra_MultiVector), intent(inout) :: this
    type(Epetra_MultiVector), intent(in) :: A
    type(error),optional,intent(out)   :: err
  end subroutine
 
  !> @name Mathematical methods
  !! @{

  !> <BR> Scales current values  this = scalar_value*this
  subroutine Scale(this,scalar_value,err)
    type(Epetra_MultiVector), intent(inout) :: this
    real(c_double),            intent(in)    :: scalar_value
    type(error),optional,intent(out)      :: err
  end subroutine
  
  !> @name Mathematical methods
  !! @{

  !> <BR> Replaces current values with scaled input  this = scalar_value*MultiVector
  subroutine Scale(this,scalar_value,MultiVector,err)
    type(Epetra_MultiVector), intent(inout) :: this
    type(Epetra_MultiVector), intent(in)    :: MultiVector 
    real(c_double),            intent(in)    :: scalar_value
    type(error),optional,intent(out)      :: err
  end subroutine


  !> @name Mathematical methods
  !! @{

  !> <BR> Computes 1-norm of each vector in the input
  function Norm1(this,err) result(Norm1_val)
    type(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: Norm1_val 
  end function 
 
  !> @name Mathematical methods
  !! @{

  !> <BR> Computes 2-norm of each vector in the input
  function Norm2(this,err) result(Norm2_val)
    type(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: Norm2_val 
  end function
 
  !> @name Mathematical methods
  !! @{

  !> <BR> Computes infinity norm of each vector in the input
  function NormInf(this,err) result(NormInf_val)
    type(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: NormInf_val 
  end function 

  !> @name Mathematical methods
  !! @{

  !> <BR> Computes weighted  norm (RMS norm) of each vector in the input
  function NormWeighted(this,weights,err) result(NormWeighted_val)
    type(Epetra_MultiVector)   ,intent(in)  :: this
    type(Epetra_MultiVector)   ,intent(in)  :: weights 
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: NormWeighted_val 
  end function 

  !> @name Mathematical methods
  !! @{

  !> <BR> Computes minimum value of each vector in the input
  function MinValue(this,err) result(MinValue_val)
    type(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: MinValue_val
  end function
 
  !> @name Mathematical methods
  !! @{

  !> <BR> MaxValue: compute maximum value of each vector in the input
 function MaxValue(this,err) result(MaxValue_val)
    type(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: MaxValue_val
  end function

  !> @name Mathematical methods
  !! @{

  !> <BR> Computes mean (average) value of each vector in the input
 function MeanValue(this,err) result(MeanValue_val)
    type(Epetra_MultiVector)   ,intent(in)  :: this
    type(error) ,optional    ,intent(out) :: err
    real(c_double) ,dimension(:),allocatable :: MeanValue_val
  end function

  !> @name Mathematical methods
  !! @{

  !> <BR> Matrix-matrix multiplication  this = ScalarThis*This + ScalarAB*MATMUL(A,B)
  subroutine Multiply(this,TransA,TransB,ScalarAB,A,B,ScalarThis,err) 
    type(Epetra_MultiVector)   ,intent(in) :: this
    character(c_char), intent(in) :: TransA
    character(c_char), intent(in) :: TransB
    type(Epetra_MultiVector)   ,intent(in)  :: A,B
    real(c_double), intent(in) :: ScalarAB, ScalarThis 
    type(error) ,optional    ,intent(out) :: err
  end subroutine

  !> @name Mathematical methods
  !! @{

  !> <BR> Element-by-element multiplication  this = ScalarThis*This + ScalarAB*A*B
  subroutine Multiply_ByEl(this,ScalarAB,A,B,ScalarThis,err)
    type(Epetra_MultiVector)   ,intent(in) :: this
    type(Epetra_MultiVector)   ,intent(in)  :: A,B
    real(c_double), intent(in) :: ScalarAB, ScalarThis
    type(error) ,optional    ,intent(out) :: err
  end subroutine

  !> @name Mathematical methods
  !! @{

  !> <BR> Element-by-element multiplication by reciprocal   this = ScalarThis*This + ScalarAB*A/B
  subroutine ReciprocalMultiply(this,ScalarAB,A,B,ScalarThis,err)
    type(Epetra_MultiVector)   ,intent(in) :: this
    type(Epetra_MultiVector)   ,intent(in)  :: A,B
    real(c_double), intent(in) :: ScalarAB, ScalarThis
    type(error) ,optional    ,intent(out) :: err
  end subroutine


  !> @name Mathematical methods
  !! @{

  !> <BR> Updates with scaled copy of input this = ScalarThis*This + ScalarA*A
  subroutine Update(this,scalarA,A,scalarThis,err)
    type(Epetra_MultiVector) ,intent(inout) :: this
    real(c_double)            ,intent(in) :: scalarA
    type(Epetra_MultiVector) ,intent(in) :: A
    real(c_double)            ,intent(in) :: scalarThis
    type(error) ,optional  ,intent(out):: err
  end subroutine 

  !> @name Mathematical methods
  !! @{

  !> <BR> Updates with scaled copies of inputs this = ScalarThis*This + ScalarA*A + ScalarB*B
  subroutine Update(this,scalarA,A,scalarB,B,scalarThis,err)
    type(Epetra_MultiVector) ,intent(inout) :: this
    real(c_double)            ,intent(in) :: scalarA
    type(Epetra_MultiVector) ,intent(in) :: A
    real(c_double)            ,intent(in) :: scalarB
    type(Epetra_MultiVector) ,intent(in) :: B
    real(c_double)            ,intent(in) :: scalarThis
    type(error) ,optional  ,intent(out):: err
  end subroutine


  !> @name Extraction methods
  !! @{

  !> <BR> Copies multivector contents into target
  function ExtractCopy(this,MyLDA,err) result(A)
    type(Epetra_MultiVector),intent(in) :: this
    real(c_double),dimension(:,:),allocatable :: A
    integer(c_int),intent(in) :: MyLDA
    type(error),optional,intent(out) :: err
  end function


  !> @name Attribute access
  !! @{

  !> <BR> Number of vectors in  multivector
  integer(c_int) function NumVectors(this)
    type(Epetra_MultiVector) ,intent(in) :: this
  end function 

  !> @name Attribute access
  !! @{

  !> <BR> Local vector length
  integer(c_int) function MyLength(this)
    type(Epetra_MultiVector) ,intent(in) :: this
  end function 

  !> @name Attribute access
  !! @{

  !> <BR> Global vector length
  integer(c_int) function GlobalLength(this)
    type(Epetra_MultiVector) ,intent(in) :: this
  end function 

  !> @name Attribute access
  !! @{

  !> <BR> Stride between successive vectors in multivector (only meaningful if ConstantStride()==true)
  integer(c_int) function Stride(this)
    type(Epetra_MultiVector) ,intent(in) :: this
  end function 

  !> @name Attribute access
  !! @{

  !> <BR> True if stride between successive vectors is constant
  integer(FT_boolean_t) function ConstantStride(this)
    use ForTrilinos_enums ,only:FT_Epetra_MultiVector_ID_t,FT_boolean_t
    type(Epetra_MultiVector) ,intent(in) :: this
  end function 

! !$ subroutine Export_UsingExporter(this,A,exporter,CombineMode,indexor,err)
! !$    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
! !$    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
! !$    use FEpetra_Export, only: Epetra_Export
! !$    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
! !$    type(Epetra_MultiVector), intent(in) :: this,A 
! !$    type(Epetra_Export),intent(in) :: exporter
! !$    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
! !$    type(Epetra_OffsetIndex), intent(in) :: indexor
! !$    type(error),optional,intent(out) :: err 
! !$    integer(c_int)     :: error_out
! !$    call this%DistObject%export(A%DistObject,exporter,CombineMode,indexor,err)
! !$    if (present(err)) err=error(error_out,'Epetra_MultiVector%Export_UsingExporter: failed')
! !$  end subroutine
! !$
! !$  subroutine Export_UsingImporter(this,A,importer,CombineMode,indexor,err)
! !$    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
! !$    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
! !$    use FEpetra_Import, only: Epetra_Import
! !$    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
! !$    type(Epetra_MultiVector), intent(in) :: this,A
! !$    type(Epetra_Import),intent(in) :: importer
! !$    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
! !$    type(Epetra_OffsetIndex), intent(in) :: indexor
! !$    type(error),optional,intent(out) :: err 
! !$    integer(c_int)     :: error_out
! !$    call this%DistObject%export(A%DistObject,importer,CombineMode,indexor,err)
! !$    if (present(err)) err=error(error_out,'Epetra_MultiVector%Export_UsingImporter: failed')
! !$  end subroutine
! !$
! !$ subroutine Import_UsingExporter(this,A,exporter,CombineMode,indexor,err)
! !$    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
! !$    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
! !$    use FEpetra_Export, only: Epetra_Export
! !$    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
! !$    type(Epetra_MultiVector), intent(in) :: this,A 
! !$    type(Epetra_Export),intent(in) :: exporter
! !$    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
! !$    type(Epetra_OffsetIndex), intent(in) :: indexor
! !$    type(error),optional,intent(out) :: err 
! !$    integer(c_int)     :: error_out
! !$    call this%DistObject%import(A%DistObject,exporter,CombineMode,indexor,err)
! !$    if (present(err)) err=error(error_out,'Epetra_MultiVector%Import_UsingExporter: failed')
! !$  end subroutine
! !$
! !$  subroutine Import_UsingImporter(this,A,importer,CombineMode,indexor,err)
! !$    use FEpetra_SrcDistObject, only: Epetra_SrcDistObject
! !$    use FEpetra_OffsetIndex, only: Epetra_OffsetIndex
! !$    use FEpetra_Import, only: Epetra_Import
! !$    use ForTrilinos_enum_wrappers, only: FT_Epetra_CombineMode_E_t
! !$    type(Epetra_MultiVector), intent(in) :: this,A
! !$    type(Epetra_Import),intent(in) :: importer
! !$    integer(FT_Epetra_CombineMode_E_t), intent(in) :: CombineMode
! !$    type(Epetra_OffsetIndex), intent(in) :: indexor
! !$    type(error),optional,intent(out) :: err 
! !$    integer(c_int)     :: error_out
! !$    call this%DistObject%import(A%DistObject,importer,CombineMode,indexor,err)
! !$    if (present(err)) err=error(error_out,'Epetra_MultiVector%Import_UsingImporter: failed')
! !$  end subroutine
! !$

end module 

