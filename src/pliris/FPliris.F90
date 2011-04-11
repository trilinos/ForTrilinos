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


module FPliris
  use ForTrilinos_enums ,only: FT_Pliris_ID_t,ForTrilinos_Universal_ID_t
  use ForTrilinos_table_man
  use ForTrilinos_universal ,only : universal
  use ForTrilinos_error ,only : error
  use iso_c_binding ,only : c_int,c_double
  use forpliris
  implicit none
  private                   ! Hide everything by default
  public :: Pliris ! Expose type/constructors/methods

  type ,extends(universal)      :: Pliris 
    private
    type(FT_Pliris_ID_t) :: FT_Pliris_id 
  contains
     !User interface -- type-bound procedures intended for use by end applications:
     procedure         :: SetLHS
     procedure         :: SetRHS
     procedure         :: SetMatrix
    !procedure         :: SetMatrix_Serial
     procedure         :: GetDistribution
     procedure         :: FactorSolve
    !procedure         :: FactorSolve_Serial
     procedure         :: Factor
     procedure         :: Solve

     !Developers only -- to be called by developers from other ForTrilinos modules, not by end applications:
     procedure         :: invalidate_id => invalidate_PlirisID
     procedure         :: ctrilinos_delete => ctrilinos_delete_Pliris
     procedure         :: get_PlirisID 
     procedure ,nopass :: alias_PlirisID
     procedure         :: generalize 
  end type

   interface Pliris ! Constructors
     !User interface -- constructors for use by end applications:
     module procedure from_scratch,from_thin_air
     !Developers only -- to be called by developers from other ForTrilinos modules, not by end applications:
     module procedure from_struct
   end interface
 
contains

  ! -------------------- User-defined constructors ------------------ 
   
  type(Pliris) function from_struct(id)
     type(FT_Pliris_ID_t) ,intent(in) :: id
     from_struct%FT_Pliris_id = id
     call from_struct%register_self
  end function

  type(Pliris) function from_thin_air()
    from_thin_air = from_struct(Pliris_Create_Default())
  end function

  type(Pliris) function from_scratch(A,x,b)
    use FEpetra_Vector ,only : Epetra_Vector
    use FEpetra_MultiVector ,only : Epetra_MultiVector
    type(Epetra_Vector) ,intent(in) :: A
    type(Epetra_MultiVector) ,intent(in) :: x,b
    type(FT_Pliris_ID_t) :: from_scratch_id
    from_scratch_id = Pliris_Create(A%get_EpetraVector_ID(),x%get_EpetraMultiVector_ID(),b%get_EpetraMultiVector_ID())
    from_scratch = from_struct(from_scratch_id)
  end function

  !----------------- Destructor ---------------------------------------

  subroutine ctrilinos_delete_Pliris(this)
    class(Pliris),intent(inout) :: this
    call Pliris_Destroy( this%FT_Pliris_id ) 
  end subroutine

  
  !----------------- Data access ---------------------------------------------

  type(FT_Pliris_ID_t) function get_PlirisID(this)
    class(Pliris) ,intent(in) :: this 
    get_PlirisID=this%FT_Pliris_id
  end function

  !----------------- Type casting ---------------------------------------------
  
  type(FT_Pliris_ID_t) function alias_PlirisID(generic_id)
    use ForTrilinos_table_man, only : CT_Alias
    use ForTrilinos_enums ,only: ForTrilinos_Universal_ID_t,FT_Pliris_ID
    use iso_c_binding     ,only: c_loc,c_int
    type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
    type(ForTrilinos_Universal_ID_t) ,allocatable ,target :: alias_id
    integer(c_int) :: status
    type(error) :: ierr
    allocate(alias_id,source=CT_Alias(generic_id,FT_Pliris_ID),stat=status)
    ierr=error(status,'FPliris: allocation failed in alias_PlirisID')
    call ierr%check_success()
    alias_PlirisID=degeneralize_Pliris(c_loc(alias_id))
  end function

  type(ForTrilinos_Universal_ID_t) function generalize(this)
   ! ____ Use for ForTrilinos function implementation ______
   use ForTrilinos_utils ,only: generalize_all
   use iso_c_binding     ,only: c_loc
   class(Pliris) ,intent(in) ,target :: this
   generalize = generalize_all(c_loc(this%FT_Pliris_id))
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !class(Pliris) ,intent(in) ,target :: this
   !generalize = Pliris_Generalize ( this%FT_Pliris_id)
   ! ____ Use for CTrilinos function implementation ______
  end function

 type(FT_Pliris_ID_t) function degeneralize_Pliris(generic_id) bind(C)
   ! ____ Use for ForTrilinos function implementation ______
    use ForTrilinos_enums ,only : FT_Pliris_ID_t
    use ,intrinsic :: iso_c_binding ,only: c_ptr,c_f_pointer
    type(c_ptr)                   ,value   :: generic_id
    type(FT_Pliris_ID_t) ,pointer :: local_ptr=>null()
    call c_f_pointer (generic_id, local_ptr)
    degeneralize_Pliris = local_ptr
   ! ____ Use for ForTrilinos function implementation ______

   ! ____ Use for CTrilinos function implementation ______
   !use ForTrilinos_enums ,only : ForTrilinos_Universal_ID_t,FT_Pliris_ID_t
   !type(ForTrilinos_Universal_ID_t) ,intent(in) :: generic_id
   !degeneralize_Pliris = Pliris_Degeneralize(generic_id)
   ! ____ Use for CTrilinos function implementation ______
  end function
 
  subroutine invalidate_PlirisID(this)
    class(Pliris),intent(inout) :: this
    this%FT_Pliris_id%table = FT_Invalid_ID
    this%FT_Pliris_id%index = FT_Invalid_Index 
    this%FT_Pliris_id%is_const = FT_FALSE
  end subroutine

  subroutine SetLHS(this,x)  
    use FEpetra_MultiVector ,only : Epetra_MultiVector
    class(Pliris) ,intent(in) :: this
    type(Epetra_MultiVector) ,intent(inout) :: x
    integer(c_int) :: success
    success = Pliris_SetLHS ( this%FT_Pliris_ID, x%get_EpetraMultiVector_ID() ) 
  end subroutine

  subroutine SetRHS (this,b) 
    use FEpetra_MultiVector ,only : Epetra_MultiVector
    class(Pliris) ,intent(in) :: this
    type(Epetra_MultiVector) ,intent(inout) :: b
    integer(c_int) :: success
    success = Pliris_SetRHS ( this%FT_Pliris_ID, b%get_EpetraMultiVector_ID() ) 
  end subroutine

  subroutine SetMatrix(this,A)
    use FEpetra_Vector ,only : Epetra_Vector
    class(Pliris) ,intent(in) :: this
    type(Epetra_Vector) ,intent(inout) :: A
    integer(c_int) :: success
    success = Pliris_SetMatrix(this%FT_Pliris_ID,A%get_EpetraVector_ID()) 
  end subroutine 

! subroutine SetMatrix_Serial (this,A)
!   use FEpetra_SerialDenseVector ,only : Epetra_SerialDenseVector
!   class(Pliris) ,intent(in) :: this
!   type(Epetra_SerialDenseVector) ,intent(inout) :: A
!   integer(c_int) :: success
!   success = Pliris_SetMatrix_Serial(this%FT_Pliris_ID,A%get_EpetraSerialDenseVector_ID()) 
! end subroutine 

  subroutine GetDistribution ( this, nprocs_row, number_of_unknowns, nrhs, my_rows, &
        my_cols, my_first_row, my_first_col, my_rhs, my_row, my_col ) 
    class(Pliris)  ,intent(in)  :: this
    integer(c_int) ,intent(in)  :: nprocs_row,number_of_unknowns,nrhs
    integer(c_int) ,intent(out) :: my_rows,my_cols,my_first_row,my_first_col,my_rhs,my_row,my_col
    integer(c_int)              :: success
    success = Pliris_GetDistribution( this%FT_Pliris_id, nprocs_row, number_of_unknowns, nrhs, my_rows, &
        my_cols, my_first_row, my_first_col, my_rhs, my_row, my_col )
  end subroutine

  subroutine FactorSolve ( this, A, my_rows, my_cols, matrix_size, num_procsr, num_rhs, secs )
    use FEpetra_Vector ,only : Epetra_Vector
    class(Pliris)       ,intent(in)  :: this
    type(Epetra_Vector) ,intent(in)  ,value :: A
    integer(c_int)      ,intent(in)  ,value :: my_rows,my_cols
    integer(c_int)      ,intent(in)  :: matrix_size,num_procsr,num_rhs
    real(c_double)      ,intent(out) :: secs
    integer(c_int) :: success
    success = Pliris_FactorSolve(this%FT_Pliris_id,A%get_EpetraVector_ID(),my_rows,my_cols,matrix_size,num_procsr,num_rhs,secs)
  end subroutine

! function Pliris_FactorSolve_Serial ( this, A, my_rows, my_cols, matrix_size, num_procsr, num_rhs, secs )
!   class(Pliris)                  ,intent(in)  :: this
!   type(Epetra_SerialDenseVector) ,intent(in)  :: A
!   integer(c_int)                 ,intent(in)  :: my_rows,my_cols,matrix_size,num_procsr,num_rhs
!   real(c_double)                 ,intent(out) :: secs
!   success = Pliris_FactorSolve_Serial( &
!     this%FT_Pliris_id,A%get_EpetraSerialDenseVector_ID(),my_rows,my_cols,matrix_size,num_procsr,num_rhs,secs)
! end function

  subroutine Factor ( this, A, matrix_size, num_procsr, permute, secs ) 
    use FEpetra_Vector ,only : Epetra_Vector
    class(Pliris)       ,intent(in)  :: this
    type(Epetra_Vector) ,intent(in)  :: A
    integer(c_int)      ,intent(in)  :: matrix_size,num_procsr
    integer(c_int)      ,intent(out) ,dimension(*) :: permute
    real(c_double)      ,intent(out)               :: secs
    integer(c_int) :: success
    success = Pliris_Factor ( this%FT_Pliris_id, A%get_EpetraVector_id(), matrix_size, num_procsr, permute, secs ) 
  end subroutine

  subroutine Solve ( this, permute, num_rhs ) 
    class(Pliris)  ,intent(in)   :: this
    integer(c_int) ,intent(in)   ,dimension(*) :: permute
    integer(c_int) ,intent(in)   :: num_rhs
    integer(c_int) :: success
    success = Pliris_Solve ( this%FT_Pliris_id, permute, num_rhs )
  end subroutine
end module 
