module ForTrilinos_hermetic_new
  private
  public :: hermetic
  type hermetic
      private
      integer, pointer :: count => null()
  contains
      procedure, private, non_overridable :: grab
      procedure, private, non_overridable :: release
      procedure, private :: assign
      procedure, private :: remote_dealloc
      procedure, non_overridable :: force_finalize
      final :: finalize_hermetic
      generic :: assignment(=) => assign
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

  subroutine assign (lhs, rhs)
      class (hermetic), intent(inout) :: lhs
      class (hermetic), intent(in) :: rhs
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
      allocate (constructor%count, source=0)
      call constructor%grab
  end function
end module
