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

module ForTrilinos_error
  use ForTrilinos_enums
  use iso_c_binding, only:c_int,c_char,c_double
  use ForTrilinos_assertion_utility ,only : error_message,assert,assert_identical

  implicit none
  private
  public :: error
 !public :: deallocate_and_check_error ! As of 3/3/2011, this is not used anywhere in ForTrilinos: we have no explicit deallocations
                                       ! because we allocate memory via allocatable entities only (never via pointers).
  type ,extends(error_message) :: error
    private
    integer(c_int) :: code
  contains
    procedure :: error_code
    procedure :: check_success
  end type

  interface error ! constructor
    module procedure new_error  
  end interface

  interface deallocate_and_check_error
    module procedure deallocate_real_rank1,deallocate_real_rank2
    module procedure deallocate_integer_rank1,deallocate_integer_rank2
  end interface
  
contains
  function new_error(new_code,new_message) 
    integer(c_int) ,intent(in) :: new_code
    character(len=*) ,intent(in) ,optional :: new_message
    type(error) :: new_error
    new_error%code = new_code 
    if (present(new_message)) new_error%error_message = error_message(new_message)
  end function

  integer(c_int) function error_code(this)
    class(error) ,intent(in) :: this
    error_code = this%code
  end function

  subroutine check_success(this)
    class(error), intent(in) :: this
    call assert( [this%code==0], [this%error_message] )
  end subroutine

  subroutine deallocate_integer_rank1(garbage,message)
    integer(c_int), dimension(:), allocatable ,intent(out) :: garbage
    character(len=*) ,intent(in) :: message
    integer(c_int) :: status
    type(error)  :: ierr
    if (allocated(garbage)) then
      deallocate(garbage,stat=status)
      ierr=error(status,message)
      call ierr%check_success()
    endif
  end subroutine

  subroutine deallocate_integer_rank2(garbage,message)
    integer(c_int), dimension(:,:),allocatable ,intent(out) :: garbage
    character(len=*) ,intent(in) :: message
    integer(c_int) :: status
    type(error)  :: ierr
    if (allocated(garbage)) then
      deallocate(garbage,stat=status)
      ierr=error(status,message)
      call ierr%check_success()
    endif
  end subroutine

  subroutine deallocate_real_rank1(garbage,message)
    real(c_double), dimension(:),allocatable ,intent(out) :: garbage
    character(len=*) ,intent(in) :: message
    integer(c_int) :: status
    type(error) :: ierr
    if (allocated(garbage)) then
      deallocate(garbage,stat=status)
      ierr=error(status,message)
      call ierr%check_success()
    endif
  end subroutine

  subroutine deallocate_real_rank2(garbage,message)
    real(c_double), dimension(:,:) ,allocatable,intent(out) :: garbage
    character(len=*) ,intent(in) :: message
    integer(c_int) :: status
    type(error) :: ierr
    if (allocated(garbage)) then
      deallocate(garbage,stat=status)
      ierr=error(status,message)
      call ierr%check_success()
    endif
  end subroutine

end module ForTrilinos_error 

