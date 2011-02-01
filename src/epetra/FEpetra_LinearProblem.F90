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


module FEpetra_LinearProblem
  use ForTrilinos_enums !,only: FT_Epetra_LinearProblem_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_universal, only: universal
  use ForTrilinos_table_man
  use ForTrilinos_hermetic,only:hermetic
  use ForTrilinos_error
  use FEpetra_MultiVector, only: Epetra_MultiVector
  use FEpetra_RowMatrix, only: Epetra_RowMatrix
  !use FEpetra_Operator, only: Epetra_Operator
  use iso_c_binding      ,only: c_int
  use forepetra
  implicit none
  private                     ! Hide everything by default
  public :: Epetra_LinearProblem ! Expose type/constructors/methods

  type, extends(universal)      :: Epetra_LinearProblem 
    private
    type(FT_Epetra_LinearProblem_ID_t) :: LinearProblem_id
  contains
     !Developers only
     procedure         :: invalidate_id => invalidate_EpetraLinearProblem_ID
     procedure         :: ctrilinos_delete => ctrilinos_delete_EpetraLinearProblem
     procedure         :: get_EpetraLinearProblem_ID 
     procedure ,nopass :: alias_EpetraLinearProblem_ID
     procedure         :: generalize 
     !Integrity check method
     procedure         :: CheckInput
     !Computational methods
     !procedure         :: LeftScale
     !procedure         :: RightScale
     !Accessor methods
     !procedure         :: GetOperator
     !precedure         :: GetMatrix
     !procedure         :: GetLHS
     !procedure         :: GetRHS
     !procedure         :: GetPDL
     !procedure         :: IsOperatorSymmetric
     !Set methods
     procedure          :: AssertSymmetric 
     procedure         :: SetPDL
     !procedure         :: SetOperator_Operator
     procedure         :: SetOperator_Matrix
     generic    :: SetOperator=>SetOperator_Matrix!,SetOperator_Operator
     procedure         :: SetLHS
     procedure         :: SetRHS
  end type

   interface Epetra_LinearProblem ! constructors
     module procedure from_scratch,duplicate,from_struct,from_scratch_matrix!,from_scratch_operator
   end interface
 
contains

  type(Epetra_LinearProblem) function from_struct(id)
    type(FT_Epetra_LinearProblem_ID_t),intent(in) :: id
    from_struct%LinearProblem_id = id
    call from_struct%register_self
  end function

!> <BR> Original C++ prototype:
!! Epetra_LinearProblem(void);
!> <BR> <BR> CTrilinos prototype:
!! CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create (  );

  type(Epetra_LinearProblem) function from_scratch()
    type(FT_Epetra_LinearProblem_ID_t) :: from_scratch_id
    from_scratch_id = Epetra_LinearProblem_Create()
    from_scratch = from_struct(from_scratch_id)
  end function

!> <BR> Original C++ prototype:
!! Epetra_LinearProblem(Epetra_RowMatrix * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
!> <BR> <BR> CTrilinos prototype:
!! CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromMatrix ( CT_Epetra_RowMatrix_ID_t AID, CT_Epetra_MultiVector_ID_t XID, 
!CT_Epetra_MultiVector_ID_t BID );

  type(Epetra_LinearProblem) function from_scratch_matrix(A,X,B)
    class(Epetra_RowMatrix), intent(in) :: A
    class(Epetra_MultiVector), intent(in) :: X, B 
    type(FT_Epetra_LinearProblem_ID_t) :: from_scratch_matrix_id
    from_scratch_matrix_id = Epetra_LinearProblem_Create_FromMatrix(A%get_EpetraRowMatrix_ID(),X%get_EpetraMultiVector_ID(),B%get_EpetraMultiVector_ID())
    from_scratch_matrix = from_struct(from_scratch_matrix_id)
  end function

!> <BR> Original C++ prototype:
!! Epetra_LinearProblem(Epetra_Operator * A, Epetra_MultiVector * X, Epetra_MultiVector * B);
!> <BR> <BR> CTrilinos prototype:
!! CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Create_FromOperator ( CT_Epetra_Operator_ID_t AID, CT_Epetra_MultiVector_ID_t XID,
! CT_Epetra_MultiVector_ID_t BID );


! type(Epetra_LinearProblem) function from_scratch_operator(A,X,B)
!    class(Epetra_Operator), intent(in) :: A
!    class(Epetra_MultiVector), intent(in) :: X, B
!    type(FT_Epetra_LinearProblem_ID_t) :: from_scratch_operator_id
!    from_scratch_operator_id = Epetra_LinearProblem_Create_FromOperator(A%get_EpetraOperator_ID(),X%get_EpetraMultiVector_ID(),B%get_EpetraMultiVector_ID())
!    from_scratch_operator = from_struct(from_scratch_operator_id)
!  end function

 !> <BR> Original C++ prototype:
 !! Epetra_LinearProblem(const Epetra_LinearProblem& Problem);
 !> <BR> <BR> CTrilinos prototype:
 !! CT_Epetra_LinearProblem_ID_t Epetra_LinearProblem_Duplicate ( CT_Epetra_LinearProblem_ID_t ProblemID );

  type(Epetra_LinearProblem) function duplicate(this)
    type(Epetra_LinearProblem) ,intent(in) :: this
    type(FT_Epetra_LinearProblem_ID_t) :: duplicate_id
    duplicate_id = Epetra_LinearProblem_Duplicate(this%LinearProblem_id)
    duplicate = from_struct(duplicate_id)
  end function

  integer(c_int) function CheckInput(this)
    class(Epetra_LinearProblem), intent(in) :: this
    CheckInput = Epetra_LinearProblem_CheckInput(this%LinearProblem_id)
  end function
 
  subroutine AssertSymmetric(this)
    class(Epetra_LinearProblem), intent(in) :: this
    call Epetra_LinearProblem_AssertSymmetric(this%LinearProblem_id)
  end subroutine

  subroutine SetPDL(this,PDL)
    class(Epetra_LinearProblem), intent(in) :: this
    integer(FT_ProblemDifficultyLevel_E_t), intent(in) :: PDL
    call Epetra_LinearProblem_SetPDL(this%LinearProblem_id,PDL)
 end subroutine

 subroutine SetOperator_Matrix(this,A)
   class(Epetra_LinearProblem), intent(in) :: this
   class(Epetra_RowMatrix), intent(in) :: A
   call Epetra_LinearProblem_SetOperator_Matrix(this%LinearProblem_id,A%get_EpetraRowMatrix_ID())
 end subroutine

 subroutine SetLHS(this,X)
   class(Epetra_LinearProblem), intent(in) :: this
   class(Epetra_MultiVector), intent(in) :: X
   call Epetra_LinearProblem_SetLHS(this%LinearProblem_id,X%get_EpetraMultiVector_ID())
 end subroutine
  
 subroutine SetRHS(this,B)
   class(Epetra_LinearProblem), intent(in) :: this
   class(Epetra_MultiVector), intent(in) :: B
   call Epetra_LinearProblem_SetRHS(this%LinearProblem_id,B%get_EpetraMultiVector_ID())
 end subroutine

  type(FT_Epetra_LinearProblem_ID_t) function get_EpetraLinearProblem_ID(this)
    class(Epetra_LinearProblem) ,intent(in) :: this 
    get_EpetraLinearProblem_ID=this%LinearProblem_id
  end function
 
  type(FT_Epetra_LinearProblem_ID_t) function alias_EpetraLinearProblem_ID(generic_id)
    use iso_c_binding        ,only: c_loc,c_int
    use ForTrilinos_enums    ,only: ForTrilinos_Universal_ID_t, FT_Epetra_LinearProblem_ID
    use ForTrilinos_table_man,only: CT_Alias
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,pointer    :: alias_id=>null()
    integer(c_int) :: status
    type(error) :: ierr
    if (.not.associated(alias_id)) then
      allocate(alias_id,source=CT_Alias(generic_id,FT_Epetra_LinearProblem_ID),stat=status)
      ierr=error(status,'FEpetra_LinearProblem:alias_EpetraLinearProblem_ID')
      call ierr%check_success()
    endif
    alias_EpetraLinearProblem_ID=degeneralize_EpetraLinearProblem(c_loc(alias_id))
    call deallocate_and_check_error(alias_id,'FEpetra_LinearProblem:alias_EpetraLinearProblem_ID')
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all 
   use iso_c_binding     ,only: c_loc
   class(Epetra_LinearProblem) ,intent(in) ,target :: this
   generalize =generalize_all(c_loc(this%LinearProblem_id))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(Epetra_LinearProblem) ,intent(in) ,target :: this
   !generalize = Epetra_LinearProblem_Generalize ( this%LinearProblem_id )
   ! ____ Use for CTrilinos function implementation ______
  
  end function
  
  type(FT_Epetra_LinearProblem_ID_t) function degeneralize_EpetraLinearProblem(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Epetra_LinearProblem_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)              ,value   :: generic_id
    type(FT_Epetra_LinearProblem_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_EpetraLinearProblem = local_ptr
   ! ____ Use for ForTrilinos function implementation ______
   
   ! ____ Use for CTrilinos function implementation ______
   ! type(ForTrilinos_Universal_ID_t) ,intent(in) : generic_id
   ! degeneralize_EpetraLinearProblem = Epetra_LinearProblem_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

  subroutine invalidate_EpetraLinearProblem_ID(this)
    class(Epetra_LinearProblem),intent(inout) :: this
    this%LinearProblem_id%table = FT_Invalid_ID
    this%LinearProblem_id%index = FT_Invalid_Index 
    this%LinearProblem_id%is_const = FT_FALSE
  end subroutine

  subroutine ctrilinos_delete_EpetraLinearProblem(this)
    class(Epetra_LinearProblem),intent(inout) :: this
    call Epetra_LinearProblem_Destroy( this%LinearProblem_id ) 
  end subroutine

end module 

