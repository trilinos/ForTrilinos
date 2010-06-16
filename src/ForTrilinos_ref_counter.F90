module ForTrilinos_ref_counter
  use ForTrilinos_hermetic, only : hermetic
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

contains

  subroutine grab(this)
      class(ref_counter), intent(inout) :: this
      if (associated(this%count)) then
          this%count = this%count + 1
      else
!          stop 'Error in grab: count not associated'
      end if
  end subroutine

  subroutine release(this)
      class (ref_counter), intent(inout) :: this
      if (associated(this%count)) then
          this%count = this%count - 1

          if (this%count == 0) then
              call this%obj%ctrilinos_delete
              deallocate (this%count, this%obj)
          end if
      else
 !         stop 'Error in release: count not associated'
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
      call this%release
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
