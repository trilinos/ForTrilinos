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

module FEpetra_CrsMatrix
  use ForTrilinos_enums !,only : FT_Epetra_RowMatrix_ID,FT_Epetra_CrsMatrix_ID_t,ForTrilinos_Universal_ID_t,FT_boolean_t,FT_FALSE,FT_TRUE
  use ForTrilinos_hermetic,only:hermetic
  use ForTrilinos_enum_wrappers
  use ForTrilinos_table_man
  use ForTrilinos_error
  use FEpetra_RowMatrix ,only : Epetra_RowMatrix
  use FEpetra_Map     !  ,only : Epetra_Map
  use iso_c_binding     ,only : c_int,c_long,c_double,c_char
  use forepetra
  implicit none
  private                         ! Hide everything by default
  public :: Epetra_CrsMatrix ! Expose type/constructors/methods

  type ,extends(Epetra_RowMatrix)  :: Epetra_CrsMatrix !"shell"
    private
    type(FT_Epetra_CrsMatrix_ID_t) :: CrsMatrix_id 
  contains
     !Constructor
     !Developers only
     procedure         :: ctrilinos_delete
     procedure         :: get_EpetraCrsMatrix_ID 
     procedure ,nopass :: alias_EpetraCrsMatrix_ID
     procedure         :: generalize 
    !Insertion/Replace/SumInto methods
     procedure         :: InsertGlobalValues
     procedure         :: ReplaceGlobalValues
    !Transformation methods
     procedure         :: FillComplete
     !Matrix data extraction routines
     procedure         :: ExtractGlobalRowCopy
     procedure         :: NumMyRowEntries
     procedure         :: MaxNumEntries
    !procedure         :: ExtractMyRowCopy
    !procedure         :: ExtractDiagonalCopy
    !Computational Methods
     procedure         :: Multiply_Vector
     procedure         :: Multiply => Multiply_MultiVector
     !Attribute Accessor Methods
     procedure         :: RowMatrixRowMap 
     procedure         :: RowMap
     procedure         :: NumGlobalEntries
     !procedure         :: Comm
     !Local/Global ID method
     procedure         :: MyGlobalRow
  end type

   interface Epetra_CrsMatrix ! constructors
     module procedure from_scratch,duplicate,from_struct
   end interface

contains

  type(Epetra_CrsMatrix) function from_struct(id)
    type(FT_Epetra_CrsMatrix_ID_t) ,intent(in) :: id
    from_struct%CrsMatrix_id = id
    call from_struct%set_EpetraRowMatrix_ID(from_struct%alias_EpetraRowMatrix_ID(from_struct%generalize()))
  end function
 
  ! Original C++ prototype:
  ! Epetra_CrsMatrix(Epetra_DataAccess CV, const Epetra_Map& RowMap, const int* NumEntriesPerRow,
  ! bool StaticProfile = false);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Create_VarPerRow ( CT_Epetra_DataAccess_E_t CV,
  !CT_Epetra_Map_ID_t RowMapID, const int * NumEntriesPerRow, boolean StaticProfile);

  type(Epetra_CrsMatrix) function from_scratch(CV,Row_Map,NumEntriesPerRow,StaticProfile)
    use ForTrilinos_enums,         only: FT_boolean_t,FT_FALSE,FT_TRUE
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    class(Epetra_Map),              intent(in) :: Row_Map
    integer(c_int), dimension(:),   intent(in) :: NumEntriesPerRow                   
    integer(FT_boolean_t), optional            :: StaticProfile
    type(FT_Epetra_CrsMatrix_ID_t) :: from_scratch_id
    if (present(StaticProfile)) then
      from_scratch_id = Epetra_CrsMatrix_Create_VarPerRow(CV,Row_Map%get_EpetraMap_ID(),NumEntriesPerRow,StaticProfile) 
    else
      from_scratch_id = Epetra_CrsMatrix_Create_VarPerRow(CV,Row_Map%get_EpetraMap_ID(),NumEntriesPerRow,FT_FALSE)
    endif
      from_scratch = from_struct(from_scratch_id)
  end function
  
  ! Original C++ prototype:
  ! Epetra_CrsMatrix(const Epetra_CrsMatrix& RowMatrix);
  ! CTrilinos prototype:
  ! CT_Epetra_CrsMatrix_ID_t Epetra_CrsMatrix_Duplicate ( CT_Epetra_BasicRowMatrix_ID_t RowMatrixID );

  type(Epetra_CrsMatrix) function duplicate(this)
    type(Epetra_CrsMatrix) ,intent(in) :: this 
    type(FT_Epetra_CrsMatrix_ID_t) :: duplicate_id
    duplicate_id = Epetra_CrsMatrix_Duplicate(this%CrsMatrix_id)
    duplicate = from_struct(duplicate_id)
  end function

  type(FT_Epetra_CrsMatrix_ID_t) function get_EpetraCrsMatrix_ID(this)
   class(Epetra_CrsMatrix) ,intent(in) :: this 
   get_EpetraCrsMatrix_ID=this%CrsMatrix_id
  end function
  
  type(FT_Epetra_CrsMatrix_ID_t) function alias_EpetraCrsMatrix_ID(generic_id)
    use iso_c_binding        ,only: c_loc
    use ForTrilinos_enums    ,only: FT_Epetra_CrsMatrix_ID,ForTrilinos_Universal_ID_t
    use ForTrilinos_table_man,only: CT_Alias
    type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(Fortrilinos_Universal_ID_t) ,pointer    :: alias_id
    allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_CrsMatrix_ID))
    alias_EpetraCrsMatrix_ID=degeneralize_EpetraCrsMatrix(c_loc(alias_id))
    deallocate(alias_id)
  end function


  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding ,only : c_loc
   class(Epetra_CrsMatrix) ,intent(in) ,target :: this
   generalize = generalize_all( c_loc(this%CrsMatrix_id) )
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! class(Epetra_CrsMatrix) ,intent(in) ,target :: this
   ! generalize = Epetra_CrsMatrix_Generalize ( this%CrsMatrix_id ) 
   ! ____ Use for CTrilinos function implementation ______
  end function
 
 type(FT_Epetra_CrsMatrix_ID_t) function degeneralize_EpetraCrsMatrix(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_CrsMatrix_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                     ,value   :: generic_id
    type(FT_Epetra_CrsMatrix_ID_t) ,pointer :: local_ptr
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraCrsMatrix = local_ptr
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! type(Fortrilinos_Universal_ID_t) ,intent(in) :: generic_id
   ! degeneralize_EpetraCrsMatrix = Epetra_CrsMatrix_Degeneralize( generic_id )
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  subroutine InsertGlobalValues(this,GlobalRow,NumEntries,values,indices,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: GlobalRow
   integer(c_int),          intent(in) :: NumEntries
   real(c_double),dimension(:),intent(in):: values 
   integer(c_int),dimension(:),intent(in):: indices 
   type(error), optional, intent(out) :: err
   integer(c_int)                          :: error_out
   error_out=Epetra_CrsMatrix_InsertGlobalValues(this%CrsMatrix_id,GlobalRow,NumEntries,values,indices)
   if (present(err)) err=error(error_out)
  end subroutine

  subroutine ReplaceGlobalValues(this,GlobalRow,NumEntries,values,indices,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: GlobalRow
   integer(c_int),          intent(in) :: NumEntries
   real(c_double), dimension(:)        :: values 
   integer(c_int),    dimension(:)        :: indices 
   type(error), optional, intent(out) :: err
   integer(c_int)                          :: error_out
   error_out=Epetra_CrsMatrix_ReplaceGlobalValues(this%CrsMatrix_id,GlobalRow,NumEntries,values,indices)
   if (present(err)) err=error(error_out)
  end subroutine

  subroutine FillComplete(this,OptimizeDataStorage)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   class(Epetra_CrsMatrix), intent(in) :: this
   logical,optional,        intent(in) :: OptimizeDataStorage
   integer(c_int)                      :: junk
   integer(FT_boolean_t)               :: OptimizeDataStorage_in
   if (present(OptimizeDataStorage)) then
     if (OptimizeDataStorage) OptimizeDataStorage_in=FT_TRUE
     if (.not.OptimizeDataStorage) OptimizeDataStorage_in=FT_FALSE
   else
     OptimizeDataStorage_in=FT_TRUE
   endif
   junk=Epetra_CrsMatrix_FillComplete(this%CrsMatrix_id,OptimizeDataStorage_in)
  end subroutine

 subroutine ExtractGlobalRowCopy(this,GlobalRow,length,NumEntries,values,indices,err)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: GlobalRow
   integer(c_int),          intent(in) :: length 
   integer(c_int),          intent(inout) :: NumEntries
   real(c_double), dimension(:),intent(inout):: values 
   integer(c_int), dimension(:),intent(inout):: indices 
   type(error), optional, intent(out) :: err
   integer(c_int)                          :: error_out
   error_out=Epetra_CrsMatrix_ExtractGlobalRowCopy_WithIndices(this%CrsMatrix_id,GlobalRow,length,NumEntries,values,indices)
   if (present(err)) err=error(error_out)
 end subroutine

 integer(c_int) function NumMyRowEntries(this,MyRow)
   use iso_c_binding, only : c_int
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int),          intent(in) :: MyRow
   NumMyRowEntries=Epetra_CrsMatrix_NumMyEntries(this%CrsMatrix_id,MyRow)
 end function

 integer(c_int) function MaxNumEntries(this)
   use iso_c_binding, only: c_int
   class(Epetra_CrsMatrix), intent(in) :: this
   MaxNumEntries=Epetra_CrsMatrix_MaxNumEntries(this%CrsMatrix_id)
 end function
 
 subroutine Multiply_Vector(this,TransA,x,y,err)
  use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
  use FEpetra_Vector, only:Epetra_Vector
  class(Epetra_CrsMatrix), intent(in) :: this
  logical, intent(in) :: TransA
  integer(FT_boolean_t)  :: TransA_in
  class(Epetra_Vector), intent(in) :: x
  class(Epetra_Vector), intent(in) :: y 
  type(error), optional, intent(inout) :: err
  integer(c_int)                       :: error_out
  if (TransA) then
    TransA_in=FT_TRUE
  else
    TransA_in=FT_FALSE
  endif
  error_out=Epetra_CrsMatrix_Multiply_Vector(this%CrsMatrix_id,TransA_in,x%get_EpetraVector_ID(),y%get_EpetraVector_ID())    
  if (present(err)) err=error(error_out)
 end subroutine

 subroutine Multiply_MultiVector(this,TransA,x,y,err)
  use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
  use FEpetra_MultiVector, only:Epetra_MultiVector
  class(Epetra_CrsMatrix), intent(in) :: this
  logical, intent(in) :: TransA
  integer(FT_boolean_t)  :: TransA_in
  class(Epetra_MultiVector), intent(in) :: x
  class(Epetra_MultiVector), intent(in) :: y 
  type(error), optional, intent(inout) :: err
  integer(c_int)                       :: error_out
  if (TransA) then
    TransA_in=FT_TRUE
  else
    TransA_in=FT_FALSE
  endif
  error_out=Epetra_CrsMatrix_Multiply_MultiVector(this%CrsMatrix_id,TransA_in,x%get_EpetraMultiVector_ID(),y%get_EpetraMultiVector_ID())    
  if (present(err)) err=error(error_out)
 end subroutine

 logical function MyGlobalRow(this,GID)
   use ForTrilinos_enums, only: FT_boolean_t,FT_FALSE,FT_TRUE
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int), intent(in) ::GID
   integer(FT_boolean_t) :: MyGlobalRow_out
   MyGlobalRow_out=Epetra_CrsMatrix_MyGlobalRow(this%CrsMatrix_id,GID)
   if (MyGlobalRow_out==FT_TRUE) then
     MyGLobalRow=.true.
   else
     MyGLobalRow=.false.
   endif
 end function
 
 type(Epetra_Map) function RowMatrixRowMap(this)
  use ForTrilinos_enums, only:FT_Epetra_Map_ID_t
  class(Epetra_CrsMatrix), intent(in) :: this
  type(FT_Epetra_Map_ID_t)            :: RowMap_ID
  RowMap_ID=Epetra_CrsMatrix_RowMatrixRowMap(this%CrsMatrix_id) 
  RowMatrixRowMap=Epetra_Map(RowMap_ID)
 end function

 type(Epetra_Map) function RowMap(this)
  use ForTrilinos_enums, only:FT_Epetra_Map_ID_t
  class(Epetra_CrsMatrix), intent(in) :: this
  type(FT_Epetra_Map_ID_t)            :: RowMap_ID
  RowMap_ID=Epetra_CrsMatrix_RowMap(this%CrsMatrix_id) 
  RowMap=Epetra_Map(RowMap_ID)
 end function

 integer(c_int) function NumGlobalEntries(this,row)
   class(Epetra_CrsMatrix), intent(in) :: this
   integer(c_int), intent(in) :: row
   NumGlobalEntries=Epetra_CrsMatrix_NumGlobalEntries(this%CrsMatrix_id,row)
 end function

 !type(FT_Epetra_Comm_ID_t) function Comm(this)
 !function Comm(this)
 !  use FEpetra_Comm, only:Epetra_Comm
 !  use FEpetra_MpiComm, only:Epetra_MpiComm
 !  class(Epetra_CrsMatrix), intent(in) :: this
 !  class(Epetra_MpiComm),allocatable:: Comm  
 !  class(Epetra_MpiComm):: comm_out  
 !  class(Epetra_Comm),allocatable :: comm_temp
 !  allocate(Epetra_MpiComm:: comm_temp)  
 !  comm_temp=Epetra_CrsMatrix_Comm(this%CrsMatrix_id)
 !  comm_out=comm_temp
 !  allocate(Comm,source=comm_out)
 !end function

  subroutine ctrilinos_delete(this)
    class(Epetra_CrsMatrix) ,intent(inout) :: this
    call this%ctrilinos_delete_EpetraRowMatrix()
    call Epetra_CrsMatrix_Destroy( this%CrsMatrix_id ) 
  end subroutine
end module 

