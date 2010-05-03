module ForTrilinos_assertion_utility
  use iso_fortran_env ,only : error_unit  
  implicit none
  private
  public :: error_message,assert,assert_identical

  type error_message
    private
    character(:) ,allocatable :: string
  end type

  interface error_message ! constructor
    module procedure new_message
  end interface 

contains

  type(error_message) function new_message(message)
    character(len=*), intent(in) :: message
    new_message%string = message 
  end function

  subroutine assert(assertion,text)
    logical ,dimension(:) ,intent(in) :: assertion
    type(error_message) ,dimension(:) ,intent(in) :: text
    integer :: i
    logical :: any_failures 
    call assert_identical( [size(assertion),size(text)] )
    any_failures=.false.
    do i=1,size(assertion)
      if (.not. assertion(i)) then
        any_failures=.true.
        write(error_unit,*) 'Assertion failed with message: '
        if (allocated(text(i)%string)) then
          write(error_unit,*) text(i)%string
        else
          write(error_unit,*) '(no message provided).'
        end if
      end if
    end do
    if (any_failures) stop 'Execution halted on failed assertion(s)!'
  end subroutine
  subroutine assert_identical(integers)
    integer ,dimension(:) ,intent(in) :: integers
    integer :: i
    logical :: any_mismatches
    any_mismatches = .false.
    do i=2,size(integers)
      if (integers(i) /= integers(1)) then
        any_mismatches = .true.
        write(error_unit,*) &
        'Value ',i,' does not match expected value ',integers(1)
      end if
    end do
    if (any_mismatches) stop 'Execution halted on failed assertion!'
  end subroutine
end module
