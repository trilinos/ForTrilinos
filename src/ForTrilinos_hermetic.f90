module ForTilinos_hermetic
  implicit none
  private
  public :: hermetic ! Expose type and type-bound procedures

  ! This module provides a base type that all ForTrilinos types should extend
  ! in order to inherit it memory management utilities that are useful in compilers
  ! that do not yet support Fortran 2003 final subroutines (destructors).
  ! The type-bound proceures are modeled after those published by 
  ! G. W. Stewart (2003) "Memory leaks in derived types revisited,"
  ! ACM Fortran Forum 22:3, 25-27.

  type ,abstract :: hermetic
    private
    integer, pointer :: temporary => null()  
  contains
    procedure :: SetTemp   ! Mark object as temporary
    procedure :: GuardTemp ! Increment the reference count
    procedure :: CleanTemp ! Decrement ref count & destroy when zero
    procedure(destructor) ,deferred :: final_subroutine
  end type

  abstract interface
    subroutine destructor(this)
      import :: hermetic
      class(hermetic) :: this
    end subroutine
  end interface 

contains
  subroutine SetTemp(this) ! Mark object as temporary
    class(hermetic) ,intent(inout) :: this 
    if (.not. associated (this%temporary)) allocate (this%temporary) 
    this%temporary = 1 
  end subroutine 

  subroutine GuardTemp (this) ! Increment reference count
    class (hermetic) :: this 
    integer, pointer :: t 
    if (associated (this%temporary)) then
      t => this% temporary 
      t = t + 1 
    end if 
  end subroutine

  subroutine CleanTemp(this) ! Decrement reference count and free temporary object
    class(hermetic) :: this 
    integer, pointer :: t 
    if (associated(this%temporary)) then  ! If temporary, 
      t=>this%temporary 
      if (t > 1) then
        t = t - 1                         ! Then decrement reference count
      else if (t == 1) then               ! If no references left,
        call this%final_subroutine()      ! Then call destructor.
        deallocate(t) 
      end if 
    end if 
  end subroutine 
end module
