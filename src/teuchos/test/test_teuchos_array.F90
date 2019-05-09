!Copyright 2017-2018, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TeuchosArray
#include "ForTrilinosTeuchos_config.hpp"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding

  use forteuchos

  implicit none
  character(len=26), parameter :: FILENAME='test_teuchos_array.f90'
  integer, parameter :: dp = kind(0.d0)
  integer, parameter :: long = selected_int_kind(18)

  SETUP_TEST()

  ADD_SUBTEST_AND_RUN(TeuchosArray_Basic)

  ADD_SUBTEST_AND_RUN(TeuchosArrayDbl_view)
  ADD_SUBTEST_AND_RUN(TeuchosArrayInt_view)
  ADD_SUBTEST_AND_RUN(TeuchosArrayLong_view)

  TEARDOWN_TEST()
contains

  FORTRILINOS_UNIT_TEST(TeuchosArray_Basic)

    type(TeuchosArrayInt) :: arr
    integer(C_INT), parameter, dimension(6) :: test_int = (/ -1, 1, 3, 3, 5, 7 /)
    integer(C_INT), pointer, dimension(:) :: iptr => NULL()

    OUT0('Starting TeuchosArray_Basic!')

    arr = TeuchosArrayInt(test_int)
    TEST_ASSERT(c_associated(arr%swigdata%cptr))

    ! Get a view to the dataa
    iptr => arr%view()
    TEST_ASSERT(associated(iptr))
    TEST_ASSERT(size(iptr) == size(test_int))
    TEST_ASSERT(all(iptr == test_int))

    call arr%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TeuchosArray_Basic)

  FORTRILINOS_UNIT_TEST(TeuchosArrayDbl_view)
    type(TeuchosArrayDbl) :: arr
    real(dp), parameter, dimension(6) :: test_dbl = (/ -1.0_dp, 0.7_dp, 3.1_dp, 2.8_dp, 5.1_dp, 6.9_dp /)
    real(dp), pointer, dimension(:) :: dptr => NULL()
    OUT0("Starting TeuchosArrayDbl_view!")

    arr = TeuchosArrayDbl(test_dbl)
    dptr => arr%view()
    TEST_ASSERT(associated(dptr))
    TEST_ASSERT(size(dptr) == size(test_dbl))
    TEST_FLOATING_ARRAY_EQUALITY(dptr, test_dbl, epsilon(dptr(1)))

    call arr%release(); TEST_IERR()

    OUT0("Finished TeuchosArrayDbl_view!")

  END_FORTRILINOS_UNIT_TEST(TeuchosArrayDbl_view)

  FORTRILINOS_UNIT_TEST(TeuchosArrayInt_view)
    type(TeuchosArrayInt) :: arr
    integer, parameter, dimension(6) :: test_int = (/ -1, 0, 3, 3, 5, 7 /)
    integer, pointer, dimension(:) :: iptr => NULL()
    OUT0("Starting TeuchosArrayInt_view!")

    arr = TeuchosArrayInt(test_int)
    iptr => arr%view()
    TEST_ASSERT(associated(iptr))
    TEST_ASSERT(size(iptr) == size(test_int))
    TEST_ARRAY_EQUALITY(iptr, test_int)

    call arr%release(); TEST_IERR()

    OUT0("Finished TeuchosArrayInt_view!")

  END_FORTRILINOS_UNIT_TEST(TeuchosArrayInt_view)

  FORTRILINOS_UNIT_TEST(TeuchosArrayLong_view)
    type(TeuchosArrayLongLong) :: arr
    integer(long), parameter, dimension(6) :: test_long = (/ -1_long, 0_long, 3_long, 3_long, 5_long, 7_long /)
    integer(long), pointer, dimension(:) :: lptr => NULL()
    OUT0("Starting TeuchosArrayLong_view!")

    arr = TeuchosArrayLongLong(test_long)
    lptr => arr%view()
    TEST_ASSERT(associated(lptr))
    TEST_ASSERT(size(lptr) == size(test_long))
    TEST_ARRAY_EQUALITY(lptr, test_long)

    call arr%release(); TEST_IERR()

    OUT0("Finished TeuchosArrayLong_view!")

  END_FORTRILINOS_UNIT_TEST(TeuchosArrayLong_view)

end program test_TeuchosArray
