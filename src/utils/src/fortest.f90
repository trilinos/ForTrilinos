! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
module fortest

! Provides procedures to be called from macros defined in FortranTestMacros.h.
! These procedures serve a twofold purpose: 1) they are clearer and easier to
! debug when defined in this module and 2) they allow the code emitted from the
! macros to remain under the 132 character limit defined by the standard.

! In some compilers, the kinds c_size_t, c_long, and c_long_long of
! ISO_C_BINDING are the same.  This causes problems when overloading functions
! to accept one of these arguments.  For instance, consider the procedure:
!
! subroutine procedure(a)
!   integer(c_size_t) :: a
!   ...
! end subroutine procedure
!
! If c_size_t and c_long are the same kind, then the following would be seen by
! the compiler as identical
!
! subroutine procedure(a)
!   integer(c_long) :: a
!   ...
! end subroutine procedure
!
! Consequently, procedure could not be overloaded for these two kinds of
! integers.  To circumvent this problem, macros <KIND>_IS_NOT_SAME_AS_LONG
! protect interfaces from multiply defining functions.

#include "ForTrilinosUtils_config.hpp"
#include "DBCF.h"

use DBCF_M
use iso_c_binding
use iso_fortran_env

#ifdef HAVE_MPI
use mpi
#endif

implicit none

! Private variables for tracking progress of tests
integer, private :: local_error(1), global_error(1)
integer, private :: comm_rank, failed_tests, total_tests
logical, private :: setup_called = .false.

interface fortest_array_equality
  module procedure fortest_array_equality_double1, &
                   fortest_array_equality_double2, &
                   fortest_array_equality_int1, &
                   fortest_array_equality_int2, &
#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
                   fortest_array_equality_size_t1, &
                   fortest_array_equality_size_t2, &
#endif
#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
                   fortest_array_equality_long_long1, &
                   fortest_array_equality_long_long2, &
#endif
                   fortest_array_equality_long1, &
                   fortest_array_equality_long2
end interface

interface fortest_array_inequality
  module procedure fortest_array_inequality_double1, &
                   fortest_array_inequality_double2, &
                   fortest_array_inequality_int1, &
                   fortest_array_inequality_int2, &
#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
                   fortest_array_inequality_size_t1, &
                   fortest_array_inequality_size_t2, &
#endif
#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
                   fortest_array_inequality_long_long1, &
                   fortest_array_inequality_long_long2, &
#endif
                   fortest_array_inequality_long1, &
                   fortest_array_inequality_long2
end interface

interface fortest_equality
  module procedure fortest_equality_double1, &
                   fortest_equality_int1, &
#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
                   fortest_equality_size_t1, &
#endif
#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
                   fortest_equality_long_long1, &
#endif
                   fortest_equality_long1
end interface

interface fortest_inequality
  module procedure fortest_inequality_double1, &
                   fortest_inequality_int1, &
#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
                   fortest_inequality_size_t1, &
#endif
#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
                   fortest_inequality_long_long1, &
#endif
                   fortest_inequality_long1
end interface

interface fortest_assert
  module procedure fortest_assert1, fortest_assert2
end interface

public

contains

  ! -------------------------------------------------------------------------- !

  subroutine setup_test2()
    ! ------------------------------------------------------------------------ !
    ! Setup module variables for the test.  This should be called once at the
    ! beginning of a test program
    ! ------------------------------------------------------------------------ !
#ifdef HAVE_MPI
    integer :: ierror
#endif
    ! ------------------------------------------------------------------------ !

    if (setup_called) then
      Insist(.FALSE., "SETUP_TEST can be called only once per test program")
    end if

    ! Initialize module level variables
    local_error(1) = 0
    global_error(1) = 0
    failed_tests = 0
    total_tests = 0

    ! Initialize MPI
#ifdef HAVE_MPI
    ierror = 0
    call MPI_Init(ierror); \
    call MPI_Comm_Rank(MPI_COMM_WORLD, comm_rank, ierror)
#else
    comm_rank = 0
#endif

    setup_called = .true.

    return
  end subroutine setup_test2

  ! -------------------------------------------------------------------------- !

  subroutine teardown_test2()
    ! ------------------------------------------------------------------------ !
    ! Tears down the test program.  For MPI builds, MPI_Finalize is called.
    ! ------------------------------------------------------------------------ !
#ifdef HAVE_MPI
    integer :: ierror
#endif
    character(len=45) :: str
    ! ------------------------------------------------------------------------ !
    if (failed_tests == 0) then
      write(str, '(A,I3,A)') "100% of ", total_tests, " tests PASSED"
    else
      write(str, '(I3,A,I3,A)') failed_tests, " of ", total_tests, " tests FAILED"; \
    end if
    if (comm_rank == 0) write(0, '(A)') str
#ifdef HAVE_MPI
    call MPI_Finalize(ierror)
#endif
    if (failed_tests > 0) then
      Insist(.FALSE., str)
    end if
    return
  end subroutine teardown_test2

  ! -------------------------------------------------------------------------- !

  subroutine setup_subtest(success)
    ! ------------------------------------------------------------------------ !
    ! Setup a subtest.  Sets success = .true. and increments the total number of
    ! tests run
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    ! ------------------------------------------------------------------------ !
    success = .true.
    total_tests = total_tests + 1
  end subroutine setup_subtest

  ! -------------------------------------------------------------------------- !

  subroutine teardown_subtest(test_name, success)
    ! ------------------------------------------------------------------------ !
    ! Tears down a subtest after completion.  For MPI builds, a reduction of the
    ! error count is performed to determine global success or failure
    ! ------------------------------------------------------------------------ !
    character(len=*), intent(in) :: test_name
    logical, intent(in) :: success
#ifdef HAVE_MPI
    integer ierror
#endif
    ! ------------------------------------------------------------------------ !
    local_error(1) = merge(0, 1, success)

#ifdef HAVE_MPI
    call MPI_Allreduce(local_error, global_error, 1, &
                       MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
#else
    global_error(1) = local_error(1)
#endif

    if (global_error(1) /= 0) then
      ! Test failed on one or more processors!
      if (comm_rank == 0) write(0, '(A)') "Test "//trim(test_name)//" FAILED!"
      failed_tests = failed_tests + 1
    end if

    return

  end subroutine teardown_subtest

  ! -------------------------------------------------------------------------- !

  subroutine write_error_diagnostics(filename, lineno, signature, errstr)
    ! ------------------------------------------------------------------------ !
    ! Write error diagnostics to stderr (file unit 0)
    ! ------------------------------------------------------------------------ !
    character(len=*), intent(in) :: filename
    integer, intent(in) :: lineno
    character(len=256), intent(in) :: signature
    character(len=*), intent(in), optional :: errstr
    ! ------------------------------------------------------------------------ !
    write(0, '(A,A,A,I6)') "File '", trim(filename), "', line ", lineno
    write(0, '(A)') " " // trim(signature)
    if (present(errstr)) then
      write(0, '(A)') "Error: " // trim(errstr)
    end if
  end subroutine write_error_diagnostics

  ! -------------------------------------------------------------------------- !

  subroutine write_to_proc0(string)
    ! ------------------------------------------------------------------------ !
    ! Write string to stderr (file unit 0) on process 0
    ! ------------------------------------------------------------------------ !
    character(len=*), intent(in) :: string
    ! ------------------------------------------------------------------------ !
    if (comm_rank == 0) write(0, '(A)') trim(string)
  end subroutine write_to_proc0

  ! -------------------------------------------------------------------------- !

  subroutine fortest_ierr(success, filename, lineno, ierr, serr)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename
    integer, intent(in) :: lineno
    integer(c_int), intent(in) :: ierr
    character(len=1024), intent(in) :: serr
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    if (ierr /= 0) then
      if (comm_rank == 0) then
        write(signature, '(A)') "TEST_IERR()"
        call write_error_diagnostics(filename, lineno, signature, serr)
      end if
      success = .false.
    end if
  end subroutine fortest_ierr

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_equality_double1(success, filename, lineno, &
                                            namea, a, nameb, b, tolerance)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    real(c_double), intent(in) :: a(:), b(:)
    real(c_double), intent(in) :: tolerance
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_FLOATING_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > tolerance) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_double1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_equality_double2(success, filename, lineno, &
                                            namea, a, nameb, b, tolerance)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    real(c_double), intent(in) :: a(:), b
    real(c_double), intent(in) :: tolerance
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > tolerance) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_double2

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_equality_int1(success, filename, lineno, &
                                         namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_int), intent(in) :: a(:), b(:)
    integer(c_int), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_int1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_equality_int2(success, filename, lineno, &
                                         namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_int), intent(in) :: a(:), b
    integer(c_int), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_int2

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_equality_long1(success, filename, lineno, &
                                          namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long), intent(in) :: a(:), b(:)
    integer(c_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_long1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_equality_long2(success, filename, lineno, &
                                          namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long), intent(in) :: a(:), b
    integer(c_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_long2

  ! -------------------------------------------------------------------------- !

#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
  subroutine fortest_array_equality_long_long1(success, filename, lineno, &
                                               namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long_long), intent(in) :: a(:), b(:)
    integer(c_long_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_long_long1
#endif

  ! -------------------------------------------------------------------------- !

#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
  subroutine fortest_array_equality_long_long2(success, filename, lineno, &
                                               namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long_long), intent(in) :: a(:), b
    integer(c_long_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_long_long2
#endif

  ! -------------------------------------------------------------------------- !

#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
  subroutine fortest_array_equality_size_t1(success, filename, lineno, &
                                            namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_size_t), intent(in) :: a(:), b(:)
    integer(c_size_t), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_size_t1
#endif

  ! -------------------------------------------------------------------------- !

#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
  subroutine fortest_array_equality_size_t2(success, filename, lineno, &
                                            namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_size_t), intent(in) :: a(:), b
    integer(c_size_t), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    if (abs(maxval(a - b)) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_equality_size_t2
#endif

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_inequality_double1(success, filename, lineno, &
                                            namea, a, nameb, b, tolerance)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    real(c_double), intent(in) :: a(:), b(:)
    real(c_double), intent(in) :: tolerance
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_FLOATING_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= tolerance) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_double1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_inequality_double2(success, filename, lineno, &
                                            namea, a, nameb, b, tolerance)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    real(c_double), intent(in) :: a(:), b
    real(c_double), intent(in) :: tolerance
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= tolerance) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_double2

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_inequality_int1(success, filename, lineno, &
                                         namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_int), intent(in) :: a(:), b(:)
    integer(c_int), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_int1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_inequality_int2(success, filename, lineno, &
                                         namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_int), intent(in) :: a(:), b
    integer(c_int), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_int2

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_inequality_long1(success, filename, lineno, &
                                          namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long), intent(in) :: a(:), b(:)
    integer(c_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_long1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_array_inequality_long2(success, filename, lineno, &
                                          namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long), intent(in) :: a(:), b
    integer(c_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_long2

  ! -------------------------------------------------------------------------- !

#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
  subroutine fortest_array_inequality_long_long1(success, filename, lineno, &
                                               namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long_long), intent(in) :: a(:), b(:)
    integer(c_long_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_long_long1
#endif

  ! -------------------------------------------------------------------------- !

#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
  subroutine fortest_array_inequality_long_long2(success, filename, lineno, &
                                               namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long_long), intent(in) :: a(:), b
    integer(c_long_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_long_long2
#endif

  ! -------------------------------------------------------------------------- !

#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
  subroutine fortest_array_inequality_size_t1(success, filename, lineno, &
                                            namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_size_t), intent(in) :: a(:), b(:)
    integer(c_size_t), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_size_t1
#endif

  ! -------------------------------------------------------------------------- !

#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
  subroutine fortest_array_inequality_size_t2(success, filename, lineno, &
                                            namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_size_t), intent(in) :: a(:), b
    integer(c_size_t), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    if (abs(maxval(a - b)) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_array_inequality_size_t2
#endif

  ! -------------------------------------------------------------------------- !

  subroutine fortest_equality_double1(success, filename, lineno, &
                                      namea, a, nameb, b, tolerance)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    real(c_double), intent(in) :: a, b
    real(c_double), intent(in) :: tolerance
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_FLOATING_EQUALITY'
    if (abs(a - b) > tolerance) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_equality_double1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_equality_int1(success, filename, lineno, &
                                   namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_int), intent(in) :: a, b
    integer(c_int), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_EQUALITY'
    if (abs(a - b) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_equality_int1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_equality_long1(success, filename, lineno, &
                                    namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long), intent(in) :: a, b
    integer(c_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_EQUALITY'
    if (abs(a - b) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_equality_long1

  ! -------------------------------------------------------------------------- !

#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
  subroutine fortest_equality_long_long1(success, filename, lineno, &
                                         namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long_long), intent(in) :: a, b
    integer(c_long_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_EQUALITY'
    if (abs(a - b) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_equality_long_long1
#endif

  ! -------------------------------------------------------------------------- !

#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
  subroutine fortest_equality_size_t1(success, filename, lineno, &
                                      namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_size_t), intent(in) :: a, b
    integer(c_size_t), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_EQUALITY'
    if (abs(a - b) > zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_equality_size_t1
#endif

  ! -------------------------------------------------------------------------- !

  subroutine fortest_inequality_double1(success, filename, lineno, &
                                        namea, a, nameb, b, tolerance)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    real(c_double), intent(in) :: a, b
    real(c_double), intent(in) :: tolerance
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_FLOATING_INEQUALITY'
    if (abs(a - b) <= tolerance) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_inequality_double1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_inequality_int1(success, filename, lineno, &
                                     namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_int), intent(in) :: a, b
    integer(c_int), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_INEQUALITY'
    if (abs(a - b) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_inequality_int1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_inequality_long1(success, filename, lineno, &
                                      namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long), intent(in) :: a, b
    integer(c_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_INEQUALITY'
    if (abs(a - b) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_inequality_long1

  ! -------------------------------------------------------------------------- !

#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
  subroutine fortest_inequality_long_long1(success, filename, lineno, &
                                           namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_long_long), intent(in) :: a, b
    integer(c_long_long), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_INEQUALITY'
    if (abs(a - b) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_inequality_long_long1
#endif

  ! -------------------------------------------------------------------------- !

#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
  subroutine fortest_inequality_size_t1(success, filename, lineno, &
                                        namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    integer(c_size_t), intent(in) :: a, b
    integer(c_size_t), parameter :: zero=0
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_INEQUALITY'
    if (abs(a - b) <= zero) then
      if (comm_rank == 0) then
        signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_inequality_size_t1
#endif

  ! -------------------------------------------------------------------------- !

  subroutine fortest_assert1(success, filename, lineno, namec, c)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namec
    integer, intent(in) :: lineno
    logical, intent(in) :: c
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ASSERT'
    if (.not.(c)) then
      if (comm_rank == 0) then
        signature = name // '(' // trim(namec) // ')'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_assert1

  ! -------------------------------------------------------------------------- !

  subroutine fortest_assert2(success, filename, lineno, namec, c)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namec
    integer, intent(in) :: lineno
    logical(c_bool), intent(in) :: c
    character(len=30) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ASSERT'
    if (.not.(c)) then
      if (comm_rank == 0) then
        signature = name // '(' // trim(namec) // ')'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_assert2

  ! -------------------------------------------------------------------------- !

  subroutine fortest_throw(success, filename, lineno, namec, ierr)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namec
    integer, intent(in) :: lineno
    integer(c_int), intent(in) :: ierr
    character(len=20) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_THROW'
    if (ierr == 0) then
      if (comm_rank == 0) then
        signature = trim(name) //  '(' // trim(namec) // ')'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_throw

  ! -------------------------------------------------------------------------- !

  subroutine fortest_nothrow(success, filename, lineno, namec, ierr)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namec
    integer, intent(in) :: lineno
    integer(c_int), intent(in) :: ierr
    character(len=20) :: name
    character(len=256) :: signature
    ! ------------------------------------------------------------------------ !
    name = 'TEST_NOTHROW'
    if (ierr /= 0) then
      if (comm_rank == 0) then
        signature = trim(name) //  '(' // trim(namec) // ')'
        call write_error_diagnostics(filename, lineno, signature)
      end if
      success = .false.
    end if
    return
  end subroutine fortest_nothrow

  ! -------------------------------------------------------------------------- !

end module fortest
