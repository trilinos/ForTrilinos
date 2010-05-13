program main
#include "all_build_macros.h"
  use iso_fortran_env ,only : error_unit ,output_unit
  use TEST_CALLS_FILE(CLASS)
  implicit none

  logical :: success
  character(len=50) :: which_test
  integer :: arg_num,arg_len,stat

  success = .FALSE.
  arg_num = 1

  call get_command_argument(arg_num,which_test,arg_len,stat)
  if (stat .ne. 0) stop "Could not determine which test to run. TEST FAILED"

  print *,"Testing ",TEST_FILE_STR(CLASS),"::",trim(which_test),TEST_SUFFIX_STR()

  call test_setup()

  print *,"Running the actual test..."
  success = select_test(which_test)

  call test_teardown()

  if (success) then
    write(output_unit,*) 
    write(output_unit,fmt='(a)') "TEST PASSED" ;
  else
    write(output_unit,*) 
    write(output_unit,fmt='(a)') "TEST FAILED" ;
    stop 1
  end if

  contains

    subroutine test_setup()
      print *,"Running test_setup()..."
    end subroutine

    subroutine test_teardown()
      print *,"Running test_teardown()..."
    end subroutine
  
end program
