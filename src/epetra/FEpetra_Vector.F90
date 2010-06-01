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

module FEpetra_Vector
  use ForTrilinos_enums   ,only: FT_Epetra_MultiVector_ID_t,FT_Epetra_Vector_ID_t,FT_Epetra_BlockMap_ID_t,ForTrilinos_Universal_ID_t,FT_boolean_t 
  use ForTrilinos_table_man
  use ForTrilinos_error
  use FEpetra_MultiVector ,only: Epetra_MultiVector
  use FEpetra_BlockMap    !,only: Epetra_BlockMap !use to circumvent reported compiler bug
  use iso_c_binding       ,only: c_int
  use forepetra
  private                                    ! Hide everything by default
  public :: Epetra_Vector!,Epetra_MultiVector ! Expose type/constructors/methods
  implicit none

  type ,extends(Epetra_MultiVector)      :: Epetra_Vector !"shell"
    private
    type(FT_Epetra_Vector_ID_t)  :: vector_id 
  contains
     procedure         :: remote_dealloc
     procedure         :: get_EpetraVector_ID 
     procedure ,nopass :: alias_EpetraVector_ID
     procedure         :: generalize 
     !Post-construction modfication routines
     procedure         :: ReplaceGlobalValues
     ! Extraction methods
     procedure         :: ExtractCopy_EpetraVector
     generic :: ExtractCopy => ExtractCopy_EpetraVector
  end type

   interface Epetra_Vector ! constructors
     module procedure constructor1,constructor2,duplicate,from_struct
   end interface
 
contains

  type(Epetra_Vector) function from_struct(id)
    type(FT_Epetra_Vector_ID_t) ,intent(in) :: id
    from_struct%vector_id = id
    from_struct%Epetra_MultiVector=Epetra_MultiVector(from_struct%alias_EpetraMultiVector_ID(from_struct%generalize()))
  end function

  ! Original C++ prototype:
  ! Epetra_Vector(const Epetra_BlockMap& Map, bool zeroOut = true);
  ! CTrilinos prototype:
  ! CT_Epetra_Vector_ID_t Epetra_Vector_Create ( CT_Epetra_BlockMap_ID_t MapID, boolean zeroOut );

  type(Epetra_Vector) function constructor1(BlockMap,zero_initial)
    use ForTrilinos_enums ,only: FT_boolean_t,FT_FALSE,FT_TRUE,FT_Epetra_BlockMap_ID_t
    use FEpetra_BlockMap  ,only: Epetra_BlockMap
    class(Epetra_BlockMap) ,intent(in) :: BlockMap
    logical ,optional      ,intent(in) :: zero_initial
    integer(FT_boolean_t)              :: zero_out
    type(FT_Epetra_Vector_ID_t) :: constructor1_id
    if (present(zero_initial).and.zero_initial) then
     zero_out=FT_TRUE
    else
     zero_out=FT_FALSE
    endif
    constructor1_id = Epetra_Vector_Create(BlockMap%get_EpetraBlockMap_ID(),zero_out)
    constructor1 = from_struct(constructor1_id)
  end function
  
  type(Epetra_Vector) function constructor2(CV,BlockMap,V)
    use ForTrilinos_enum_wrappers
    use FEpetra_BlockMap  ,only: Epetra_BlockMap
    class(Epetra_BlockMap) ,intent(in) :: BlockMap
    integer(FT_Epetra_DataAccess_E_t), intent(in) :: CV
    real(c_double),dimension(:) :: V
    type(FT_Epetra_Vector_ID_t) :: constructor2_id
    constructor2_id = Epetra_Vector_Create_FromArray(CV,BlockMap%get_EpetraBlockMap_ID(),V)
    constructor2 = from_struct(constructor2_id) 
  end function

  ! Original C++ prototype:
  ! Epetra_Vector(const Epetra_Vector& Source);
  ! CTrilinos prototype:
  ! CT_Epetra_Vector_ID_t Epetra_Vector_Duplicate ( CT_Epetra_Vector_ID_t SourceID );

  type(Epetra_Vector) function duplicate(this)
    type(Epetra_Vector) ,intent(in) :: this 
    type(FT_Epetra_Vector_ID_t) :: duplicate_id
    duplicate_id = Epetra_Vector_Duplicate(this%vector_id)
    duplicate = from_struct(duplicate_id)
  end function

  type(FT_Epetra_Vector_ID_t) function get_EpetraVector_ID(this)
    class(Epetra_Vector) ,intent(in) :: this 
    get_EpetraVector_ID=this%vector_id
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

  subroutine ReplaceGlobalValues(this,NumEntries,values,indices,err)
    class(Epetra_Vector), intent(in) :: this
    integer(c_int),       intent(in) :: NumEntries
    real(c_double),dimension(:),intent(in) :: values
    integer(c_int),dimension(:),intent(in) :: indices 
    type(error),optional,intent(out) :: err
    integer(c_int)                      :: error_out
    error_out=Epetra_Vector_ReplaceGlobalValues(this%vector_id,NumEntries,values,indices)
    if (present(err)) err=error(error_out)
  end subroutine

  function ExtractCopy_EpetraVector(this,err) result(ExtractCopy_out)
   class(Epetra_Vector), intent(in) :: this
   real(c_double), dimension(:), allocatable::  ExtractCopy_out 
   type(error),optional,intent(out) :: err
   integer(c_int)                      :: error_out
   allocate(ExtractCopy_out(this%MyLength()))
   error_out = Epetra_Vector_ExtractCopy(this%vector_id,ExtractCopy_out)
   if (present(err)) err=error(error_out)
  end function 
  
  subroutine remote_dealloc(this)
    class(Epetra_Vector) ,intent(inout) :: this
    call this%remote_dealloc_EpetraMultiVector()
    call Epetra_Vector_Destroy(this%vector_id) 
  end subroutine

end module 

