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


  SETUP_TEST()

  ADD_SUBTEST_AND_RUN(TeuchosArray_Basic)

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

  END_FORTRILINOS_UNIT_TEST(TeuchosArray_Basic)

end program test_TeuchosArray
