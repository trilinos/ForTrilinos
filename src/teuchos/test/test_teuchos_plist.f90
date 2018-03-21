!Copyright 2017-2018, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TeuchosPList
#include "ForTrilinosTeuchos_config.hpp"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding

  use forteuchos

  implicit none
  character(len=26), parameter :: FILENAME='test_teuchos_plist.f90'


  SETUP_TEST()

  ADD_SUBTEST_AND_RUN(TeuchosPList_Basic)

  TEARDOWN_TEST()
contains

  FORTRILINOS_UNIT_TEST(TeuchosPList_Basic)

    type(ParameterList) :: plist, sublist, sublistref
    integer(C_INT), dimension(6) :: test_int = (/ -1, 1, 3, 3, 5, 7 /)
    real(C_DOUBLE), dimension(4) :: test_dbl = (/ 0.1d0, 1.9d0, -2.0d0, 4.0d0 /)

    integer(C_INT), allocatable :: iarr(:)
    real(C_DOUBLE), allocatable :: darr(:)
    integer(C_INT) :: ival
    real(C_DOUBLE) :: dval
    logical(C_BOOL) :: bval, true=.true. ! FIXME: can we get rid of this true somethow?
    character(kind=C_CHAR, len=:), allocatable :: sval

    OUT0('Starting TeuchosPList_Basic!')

    plist = ParameterList('myname'); TEST_IERR()
    TEST_ASSERT(c_associated(plist%swigdata%ptr))

    ! Test a function that raises an exception
    TEST_THROW(call load_from_xml(plist, 'nonexistent_path.xml'))
    TEST_ASSERT(c_associated(plist%swigdata%ptr))

    ! Get and set a vlaue
    call plist%set('myint', 4)
    ival = plist%get_integer('myint')
    TEST_EQUALITY(ival, 4)

    call plist%set('mydbl', 1.25_C_DOUBLE)
    dval = plist%get_real('mydbl')
    TEST_FLOATING_EQUALITY(dval, 1.25_C_DOUBLE, epsilon(1.0_C_DOUBLE))

    bval = .false.
    call plist%set('mybool', true)
    bval = plist%get_logical('mybool')
    TEST_ASSERT(bval)

    call plist%set('intarr', test_int)
    call plist%set('dblarr', test_dbl)

    ! FIXME: without these allocation, the test segfaults when compiled with flang
    ! According to [this](https://stackoverflow.com/questions/42140832/automatic-array-allocation-upon-assignment-in-fortran) it
    ! seems that it should be supported in F2003. So it is not clear why flang fails.
    allocate(iarr(10))
    allocate(darr(10))
    iarr = plist%get_arr_integer('intarr')
    darr = plist%get_arr_real('dblarr')

    ! Wrong parameter type
    TEST_THROW(iarr = plist%get_arr_integer('dblarr'))

    call plist%set('deleteme', 123)
    TEST_ASSERT(plist%is_parameter('deleteme'))

    call plist%remove('deleteme')
    TEST_ASSERT((.not. plist%is_parameter('deleteme')))

    TEST_ASSERT((.not. c_associated(sublist%swigdata%ptr)))
    sublist = plist%sublist('sublist')
    TEST_ASSERT(c_associated(sublist%swigdata%ptr))

    call sublist%set('anotherval', 4.0d0)

    call sublist%set('stringval', 'some string!')
    sval = sublist%get_string('stringval')
    TEST_EQUALITY('some string!', sval)

    call plist%set('sublist2', sublist)

    ! Add a sub-sublist
    sublist = sublist%sublist('subsublist')

    ! Add another parameter on the original sublist
    call sublist%set('late_arrival', 1.0d0)
    call sublist%release()

    sublistref = plist%sublist('sublist')
    call sublistref%set('added_to_copy', 1.0d0)
    OUT0('Printing copy of sublist...')
    call sublistref%print()
    call sublistref%release()

    OUT0('Printing...')
    call plist%print()

    OUT0('Saving to XML file')
    call save_to_xml(plist, 'myparams.xml')

    call plist%release()

    OUT0('Finished TeuchosPList_Basic!')

  END_FORTRILINOS_UNIT_TEST(TeuchosPList_Basic)

end program test_TeuchosPList
