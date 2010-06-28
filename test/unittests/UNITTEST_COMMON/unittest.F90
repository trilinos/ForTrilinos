program main
#include "ForTrilinos_config.h"

  use iso_fortran_env ,only : error_unit ,output_unit
  use TEST_CALLS_FILE
  implicit none
#ifdef HAVE_MPI
include 'mpif.h'
#endif

  logical :: success,fullsuccess
  character(len=100) :: which_test
  character(len=100),dimension(100) :: test_list
  integer :: arg_len,arg_cnt,stat,u
  integer :: ierr,test_num,test_cnt

#ifdef HAVE_MPI
  integer, dimension(1) :: localfail
  integer, dimension(1) :: anyfail

  call MPI_INIT(ierr)
#endif

  success = .FALSE.
  arg_cnt = command_argument_count()
  if (arg_cnt == 0) then
    write(output_unit,fmt='(a)') "No tests specified. TEST FAILED"
#ifdef HAVE_MPI
    call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    stop 1
#endif
  end if

  call get_command_argument(1,which_test,arg_len,stat)
  if (stat < 0) then
    write(output_unit,fmt='(a)') "Command line argument cropped. TEST FAILED"
#ifdef HAVE_MPI
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    stop 1
#endif
  end if
  if (stat > 0) then
    write(output_unit,fmt='(a)') "Could not retrieve first command line argument. TEST FAILED"
#ifdef HAVE_MPI
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    stop 1
#endif
  end if

  if (which_test == "-f") then
    if (arg_cnt < 2) then
      write(output_unit,fmt='(a)') "No test list file specified. TEST FAILED"
  #ifdef HAVE_MPI
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    stop 1
#endif  
    end if
    call get_command_argument(2,which_test,arg_len,stat)
    if (stat < 0) then
      write(output_unit,fmt='(a)') "Test file command line argument cropped. TEST FAILED"
#ifdef HAVE_MPI
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    stop 1
#endif
    end if
    if (stat > 0) then
      write(output_unit,fmt='(a)') "Could not retrieve command line argument for test file. TEST FAILED"
#ifdef HAVE_MPI
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    stop 1
#endif
    end if

    u = 100
    open(unit=u, iostat=ierr, file=which_test)
    if (ierr .NE. 0) then
      write(output_unit,fmt='(a)') "Error opening test list file. TEST FAILED"
#ifdef HAVE_MPI
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    stop 1
#endif
    end if
    test_cnt = 0
    do test_num = 1,100
      read(u, '(a100)') which_test
      if (which_test == "###ENDOFTESTS") exit
      test_list(test_num) = which_test
      test_cnt = test_cnt + 1
    end do
    close(unit=u)
  else
#ifdef HAVE_MPI
    write(output_unit,fmt='(a)') "You must specify a test file in MPI mode. TEST FAILED"
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    test_cnt = arg_cnt
    if (test_cnt > 100) then
      write(output_unit,fmt='(a)') "Too many tests specified. TEST FAILED"
      stop 1
    end if
    do test_num = 1,test_cnt
      call get_command_argument(test_num,which_test,arg_len,stat)
      test_list(test_num) = which_test
    end do
#endif
  end if

  fullsuccess = .TRUE.
  do test_num = 1,test_cnt
    which_test = test_list(test_num)

    print *,"Testing ",TEST_IMPLS_FILE_STR,"::",trim(which_test),"_UnitTest"

    call test_setup()

    print *,"Running the actual test..."
    success = select_test(which_test)
    call test_teardown()

#ifdef HAVE_MPI
    if (success) then
      localfail(1) = 0
    else
      localfail(1) = 1
    end if
    call MPI_AllReduce(localfail, anyfail, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierr)
    if (anyfail(1) == 0) then
      success = .TRUE.
    else
      success = .FALSE.
    end if
#endif

    fullsuccess = fullsuccess .and. success
    if (success) then
      write(output_unit,fmt='(a)') "TEST PASSED"
    else
      write(output_unit,fmt='(a)') "TEST FAILED"
#ifdef HAVE_MPI
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    stop 1
#endif
    end if
  end do

  if (arg_cnt > 1) then
    if (fullsuccess) then
      write(output_unit,fmt='(a)') "END RESULT: ALL TESTS PASSED" ;
    else
      write(output_unit,fmt='(a)') "END RESULT: SOME TESTS FAILED" ;
#ifdef HAVE_MPI
     call MPI_Abort(MPI_COMM_WORLD,1,ierr)
#else
    stop 1
#endif
    end if
  end if

#ifdef HAVE_MPI
  call MPI_FINALIZE(ierr)
#endif

  contains

    subroutine test_setup()
      print *,"Running test_setup()..."
    end subroutine

    subroutine test_teardown()
      print *,"Running test_teardown()..."
    end subroutine
  
end program
