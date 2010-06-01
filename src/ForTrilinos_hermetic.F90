module ForTrilinos_hermetic
  private
  public :: hermetic
  type, abstract :: hermetic
  contains
      procedure(free_memory), deferred :: remote_dealloc
  end type

  abstract interface
     subroutine free_memory (this)
        import :: hermetic
        class(hermetic), intent(inout) :: this
     end subroutine
  end interface
end module

