! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
module fortest

! Provides procedures to be called from macros defined in FortranTestUtilities.h.
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

#include "ForTrilinos_config.h"
use forerror

use iso_c_binding
use iso_fortran_env

#if FORTRILINOS_USE_MPI
use mpi
#endif

implicit none

! Private variables for tracking progress of tests
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
                   fortest_equality_char, &
#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
                   fortest_equality_size_t1, &
#endif
#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
                   fortest_equality_long_long1, &
#endif
                   fortest_equality_long1
end interface

interface isnan
  module procedure isnan_double1, &
                   isnan_int1, &
#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
                   isnan_size_t1, &
#endif
#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
                   isnan_long_long1, &
#endif
                   isnan_long1
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

public

contains

  ! -------------------------------------------------------------------------- !

  subroutine fortest_stop(message)
    implicit none
    character(len=*), intent(in), optional :: message
    integer, parameter :: stopcode=99329
    flush(error_unit)
    flush(output_unit)
    if (present(message)) then
      write(error_unit, '(A)') trim(message)
    end if
    stop stopcode
  end subroutine fortest_stop

  subroutine gather_success(local_success, global_success)
    implicit none
    logical, intent(in) :: local_success
    logical, intent(out) :: global_success
    integer :: local_error(1), global_error(1)
#if FORTRILINOS_USE_MPI
    integer ierror
#endif
    ! ------------------------------------------------------------------------ !

    local_error(1) = merge(0, 1, local_success)

#if FORTRILINOS_USE_MPI
    call MPI_Allreduce(local_error, global_error, 1, &
                       MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
#else
    global_error(1) = local_error(1)
#endif

    if (global_error(1) /= 0) then
      global_success = .false.
    else
      global_success = .true.
    end if

  end subroutine gather_success

  ! -------------------------------------------------------------------------- !

  subroutine setup_test2()
    ! ------------------------------------------------------------------------ !
    ! Setup module variables for the test.  This should be called once at the
    ! beginning of a test program
    ! ------------------------------------------------------------------------ !
#if FORTRILINOS_USE_MPI
    integer :: ierror
#endif
    ! ------------------------------------------------------------------------ !

    if (setup_called) then
      call fortest_stop("SETUP_TEST can be called only once per test program")
    end if

    ! Initialize module level variables
    failed_tests = 0
    total_tests = 0

    ! Initialize MPI
#if FORTRILINOS_USE_MPI
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
#if FORTRILINOS_USE_MPI
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
#if FORTRILINOS_USE_MPI
    call MPI_Finalize(ierror)
#endif
    if (failed_tests > 0) then
      call fortest_stop()
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
    integer :: local_error(1), global_error(1)
#if FORTRILINOS_USE_MPI
    integer ierror
#endif
    ! ------------------------------------------------------------------------ !
    local_error(1) = merge(0, 1, success)

#if FORTRILINOS_USE_MPI
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
    flush(error_unit)
    flush(output_unit)
    write(error_unit, '(A,I3)') "***ERROR: On Proc: ", comm_rank
    write(error_unit, '(A,A,A,I6)') &
      "*** ERROR: File: '", trim(filename), "', line ", lineno
    write(error_unit, '(A)') "*** ERROR: " // trim(signature)
    if (present(errstr)) then
      write(error_unit, '(A)') "*** ERROR: " // trim(errstr)
    end if
  end subroutine write_error_diagnostics

  ! -------------------------------------------------------------------------- !

  subroutine write_to_proc0(string)
    ! ------------------------------------------------------------------------ !
    ! Write string to stderr (file unit 0) on process 0
    ! ------------------------------------------------------------------------ !
    character(len=*), intent(in) :: string
    ! ------------------------------------------------------------------------ !
    if (comm_rank == 0) write(output_unit, '(A)') trim(string)
  end subroutine write_to_proc0

  ! -------------------------------------------------------------------------- !

  subroutine fortest_ierr(success, filename, lineno)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename
    integer, intent(in) :: lineno
    character(len=256) :: signature
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    lcl_success = (fortrilinos_ierr == 0)
    if (.not. lcl_success) then
      write(signature, '(A)') "TEST_IERR()"
      call write_error_diagnostics(filename, lineno, signature, fortrilinos_get_serr())
    end if
    call gather_success(lcl_success, success)
    return
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_FLOATING_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= tolerance)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_FLOATING_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= tolerance)
    if (any(isnan(a)) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= zero)
    if (any(isnan(a)) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= zero)
    if (any(isnan(a)) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_EQUALITY'
    lcl_success = (abs(maxval(a - b)) <= zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_FLOATING_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > tolerance)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > tolerance)
    if (any(isnan(a)) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > zero)
    if (any(isnan(a)) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > zero)
    if (any(isnan(a)) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ARRAY_INEQUALITY'
    lcl_success = (abs(maxval(a - b)) > zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_FLOATING_EQUALITY'
    lcl_success = (abs(a - b) <= tolerance)
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_EQUALITY'
    lcl_success = (abs(a - b) <= zero)
    if (isnan(a) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//')'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_EQUALITY'
    lcl_success = (abs(a - b) <= zero)
    if (isnan(a) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//')'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_EQUALITY'
    lcl_success = (abs(a - b) <= zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//')'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_EQUALITY'
    lcl_success = (abs(a - b) <= zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//')'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
    return
  end subroutine fortest_equality_size_t1
#endif

  ! -------------------------------------------------------------------------- !

  subroutine fortest_equality_char(success, filename, lineno, &
                                   namea, a, nameb, b)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namea, nameb
    integer, intent(in) :: lineno
    character(len=*), intent(in) :: a, b
    character(len=30) :: name
    character(len=256) :: signature
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_EQUALITY'
    lcl_success = (trim(a) == trim(b))
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//')'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
    return
  end subroutine fortest_equality_char

  ! -------------------------------------------------------------------------- !

  elemental logical function isnan_double1(a)
    real(c_double), intent(in) :: a
    isnan_double1 = (a /= a)
  end function isnan_double1

  ! -------------------------------------------------------------------------- !

  elemental logical function isnan_int1(a)
    integer(c_int), intent(in) :: a
    isnan_int1 = (a /= a)
  end function isnan_int1

  ! -------------------------------------------------------------------------- !

  elemental logical function isnan_long1(a)
    integer(c_long), intent(in) :: a
    isnan_long1 = (a /= a)
  end function isnan_long1

  ! -------------------------------------------------------------------------- !

#ifdef LONG_LONG_IS_NOT_SAME_AS_LONG
  elemental logical function isnan_long_long1(a)
    integer(c_long_long), intent(in) :: a
    isnan_long_long1 = (a /= a)
  end function isnan_long_long1
#endif

  ! -------------------------------------------------------------------------- !

#ifdef SIZE_T_IS_NOT_SAME_AS_LONG
  elemental logical function isnan_size_t1(a)
    integer(c_size_t), intent(in) :: a
    isnan_size_t1 = (a /= a)
  end function isnan_size_t1
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_FLOATING_INEQUALITY'
    lcl_success = (abs(a - b) > tolerance)
    if (isnan(a) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_INEQUALITY'
    lcl_success = (abs(a - b) > zero)
    if (isnan(a) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_INEQUALITY'
    lcl_success = (abs(a - b) > zero)
    if (isnan(a) .or. isnan(b)) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_INEQUALITY'
    lcl_success = (abs(a - b) > zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
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
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_INEQUALITY'
    lcl_success = (abs(a - b) > zero)
    if (any(isnan(a)) .or. any(isnan(b))) lcl_success = .false.
    if (.not. lcl_success) then
      signature = trim(name)//'('//trim(namea)//', '//trim(nameb)//', TOL)'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
    return
  end subroutine fortest_inequality_size_t1
#endif

  ! -------------------------------------------------------------------------- !

  subroutine fortest_assert(success, filename, lineno, namec, c)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namec
    integer, intent(in) :: lineno
    logical, intent(in) :: c
    character(len=30) :: name
    character(len=256) :: signature
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_ASSERT'
    lcl_success = c
    if (.not. lcl_success) then
      signature = trim(name) // '(' // trim(namec) // ')'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
    return
  end subroutine fortest_assert

  ! -------------------------------------------------------------------------- !

  subroutine fortest_throw(success, filename, lineno, namec)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namec
    integer, intent(in) :: lineno
    character(len=20) :: name
    character(len=256) :: signature
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_THROW'
    lcl_success = (fortrilinos_ierr /= 0)
    if (.not. lcl_success) then
      signature = trim(name) //  '(' // trim(namec) // ')'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
    fortrilinos_ierr = 0
    return
  end subroutine fortest_throw

  ! -------------------------------------------------------------------------- !

  subroutine fortest_nothrow(success, filename, lineno, namec)
    ! ------------------------------------------------------------------------ !
    logical, intent(inout) :: success
    character(len=*), intent(in) :: filename, namec
    integer, intent(in) :: lineno
    character(len=20) :: name
    character(len=256) :: signature
    logical :: lcl_success
    ! ------------------------------------------------------------------------ !
    name = 'TEST_NOTHROW'
    lcl_success = (fortrilinos_ierr == 0)
    if (.not. lcl_success) then
      signature = trim(name) //  '(' // trim(namec) // ')'
      call write_error_diagnostics(filename, lineno, signature)
    end if
    call gather_success(lcl_success, success)
    return
  end subroutine fortest_nothrow

  ! -------------------------------------------------------------------------- !

end module fortest
