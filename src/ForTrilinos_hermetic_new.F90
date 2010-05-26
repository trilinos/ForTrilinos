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

module ForTrilinos_hermetic_new
  private
  public :: hermetic
  type hermetic
      private
      integer, pointer :: count => null()
  contains
      procedure, private, non_overridable :: grab
      procedure, private, non_overridable :: release
      procedure, private :: assign_hermetic
      procedure, private :: remote_dealloc
      procedure, non_overridable :: force_finalize
      final :: finalize_hermetic
      generic :: assignment(=) => assign_hermetic
  end type

  interface hermetic
      module procedure constructor
  end interface

contains

  subroutine remote_dealloc (this)
    class(hermetic), intent(inout) :: this

    print *, 'this function provides a hook, you should never seen it invoked'
    stop 'Error: hermetic::remote_dealloc is invoked'
  end subroutine

  subroutine force_finalize(this)
    class(hermetic), intent(inout) :: this

    call this%release
  end subroutine

  subroutine grab(this)
      class(hermetic), intent(inout) :: this
      if (associated(this%count)) then
          this%count = this%count + 1
          print *, 'grab: increase count by 1', this%count
      else
          stop 'Error in grab: count not associated'
      end if
  end subroutine

  subroutine release(this)
      class (hermetic), intent(inout) :: this
      if (associated(this%count)) then
          this%count = this%count - 1
          print *, 'release: decrease count by 1', this%count

          if (this%count == 0) then
              print *, 'delete remote data'
              deallocate (this%count)
              call this%remote_dealloc
          end if
      else
          stop 'Error in release: count not associated'
      end if
  end subroutine
  
  subroutine assign_hermetic (lhs, rhs)
      class (hermetic), intent(inout) :: lhs
      class (hermetic), intent(in) :: rhs
      print *,'hermetic assign'
      if (associated(lhs%count)) call lhs%release
      lhs%count => rhs%count
      call lhs%grab
  end subroutine

  subroutine finalize_hermetic (this)
      type(hermetic), intent(inout) :: this
      call this%release
  end subroutine

  function constructor ()
      type(hermetic), allocatable :: constructor
      allocate (constructor)
      print *,'hermetic constructor'
      allocate (constructor%count, source=0)
      call constructor%grab
  end function
end module
