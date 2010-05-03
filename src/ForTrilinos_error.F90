module ForTrilinos_error
  use iso_c_binding, only:c_int,c_char
  implicit none
  private
  public :: error

  type error
    private
    integer(c_int) :: code
    character(:) ,allocatable :: message
  contains
    procedure :: error_message
    procedure :: error_code
  end type

  interface error ! constructor
    module procedure new_error  
  end interface
  
contains
  function new_error(new_code,new_message) 
    integer(c_int) ,intent(in) :: new_code
    character(:),allocatable ,intent(in) ,optional :: new_message
    type(error) :: new_error
    new_error%code = new_code 
    if (present(new_message)) new_error%message = new_message
  end function

  function error_message(this)
    class(error) ,intent(in) :: this 
    character(:) ,allocatable :: error_message
    if (allocated(this%message)) then
       error_message = this%message
    else
       error_message = 'No error message provided.'
    end if
  end function

  integer(c_int) function error_code(this)
    class(error) ,intent(in) :: this
    error_code = this%code
  end function
end module ForTrilinos_error 

