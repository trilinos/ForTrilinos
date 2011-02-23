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
! Questions? Contact Karla Morris  (knmorri@sandia.gov)
!                    Damian Rouson (rouson@sandia.gov)
!*********************************************************************

module ForTrilinos_ref_counter
  use ForTrilinos_hermetic ,only : hermetic
  use ForTrilinos_assertion_utility ,only : assert,error_message
  implicit none
  private
  public :: ref_counter
  type ref_counter
      private
      integer, pointer :: count => null()
      class(hermetic), pointer :: obj => null()
  contains
      procedure, non_overridable :: grab
      procedure, non_overridable :: release
      procedure :: assign
      final :: finalize_ref_counter
      generic :: assignment(=) => assign
  end type

  interface ref_counter
      module procedure constructor
  end interface

#ifdef ForTrilinos_ASSERTIONS
  logical ,parameter :: assertions=.true.
#else
  logical ,parameter :: assertions=.false.
#endif

contains

  subroutine grab(this)
    class(ref_counter), intent(inout) :: this
    if (assertions) call assert( [associated(this%count)], [error_message('Ref_counter%grab: count not associated.')] )
    this%count = this%count + 1
  end subroutine

  recursive subroutine release(this)
    use ForTrilinos_error ,only : error
    class (ref_counter), intent(inout) :: this
    integer :: status
    type(error) :: ierr
    if (assertions) call assert( [associated(this%count)], [error_message('Ref_counter%release: count not associated.')] )
    this%count = this%count - 1
    if (this%count == 0) then
      deallocate (this%count,stat=status)
      ierr=error(status,'Ref_counter%release: this%count')
      call ierr%check_success()
      deallocate (this%obj,stat=status)
      ierr=error(status,'Ref_counter%release: this%obj')
      call ierr%check_success()
    end if
  end subroutine

  subroutine assign (lhs, rhs)
    class (ref_counter), intent(inout) :: lhs
    class (ref_counter), intent(in) :: rhs
    if (associated(lhs%count)) call lhs%release
    lhs%count => rhs%count
    lhs%obj => rhs%obj
    call lhs%grab
  end subroutine

  subroutine finalize_ref_counter (this)
    type(ref_counter), intent(inout) :: this
    if (associated(this%count)) call this%release
  end subroutine

  function constructor (object)
    class(hermetic), intent(in) :: object
    type(ref_counter), allocatable :: constructor
    allocate (constructor)
    allocate (constructor%count, source=0)
    allocate (constructor%obj, source=object)
    call constructor%grab
  end function
end module
