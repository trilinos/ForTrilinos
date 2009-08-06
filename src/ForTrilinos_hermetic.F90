module ForTrilinos_hermetic
! Can't use Fortran-style include due to  presence of C pre-processor directives in ForTrilinos_config.h:
!#include "ForTrilinos_config.h"
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
    integer, pointer :: temporary => null() ! Null symbolizes "non-temporary" data
  contains
    procedure :: SetTemp      ! Mark object as temporary
    procedure :: GuardTemp    ! Increment the reference count
    procedure :: CleanTemp    ! Decrement ref count & destroy when zero
    procedure :: is_temporary ! Check for temporary mark
    procedure(destructor) ,deferred :: final_subroutine
#ifdef FINALIZATION_SUPPORTED
   final :: finalize_hermetic
#endif 
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
    if (.not. associated (this%temporary)) allocate(this%temporary) 
    this%temporary = 1 
  end subroutine 

  subroutine GuardTemp (this) 
    class (hermetic) :: this 
    integer, pointer :: t 
    if (this%is_temporary()) then ! If temporary,
      t => this% temporary 
      t = t + 1                   ! then increment count.
    end if 
  end subroutine

  subroutine CleanTemp(this) 
    class(hermetic) :: this 
    integer, pointer :: t 
    if (this%is_temporary()) then         ! If temporary, 
      t=>this%temporary 
      if (t > 1) then
        t = t - 1                         ! then decrement reference count
      else if (t == 1) then               ! else if no references left,
        call this%final_subroutine()      ! then call destructor.
        deallocate(t) 
      end if 
    end if 
  end subroutine 

  logical function is_temporary(this)
    class(hermetic) ,intent(in):: this
    if (associated(this%temporary)) then  
      is_temporary = .true.
    else
      is_temporary = .false.
    end if
  end function
end module
