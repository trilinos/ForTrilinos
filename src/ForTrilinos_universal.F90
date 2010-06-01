module ForTrilinos_universal
  use ForTrilinos_hermetic ,only : hermetic
  use ForTrilinos_ref_counter, only : ref_counter
  implicit none

  type ,abstract ,extends(hermetic) :: universal
    private
    type(ref_counter) :: counter

    contains
    procedure, non_overridable :: force_finalize
    procedure, non_overridable :: register_self
  end type

  contains

  subroutine force_finalize (this)
    class(universal), intent(inout) :: this

    call this%counter%release
  end subroutine

  subroutine register_self (this)
    class(universal), intent(inout) :: this

    this%counter = ref_counter(this)
  end subroutine
end module
