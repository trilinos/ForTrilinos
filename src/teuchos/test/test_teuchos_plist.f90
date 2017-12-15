!Copyright 2017, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TeuchosPList
#include "FortranTestUtilities.h"
#include "ForTrilinosTeuchos_config.hpp"
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
    integer, dimension(6) :: test_int = (/ -1, 1, 3, 3, 5, 7 /)
    real(C_DOUBLE), dimension(4) :: test_dbl = (/ 0.1d0, 1.9d0, -2.0d0, 4.0d0 /)

    integer :: ival
    real(C_DOUBLE) :: dval
    logical(C_BOOL) :: bval, true=.true., false=.false. ! FIXME: can we get rid of this true somethow?
    character(kind=C_CHAR, len=16) :: sval

    OUT0('Starting TeuchosPList_Basic!')

    call plist%create('myname'); TEST_IERR()
    TEST_EQUALITY_CONST(c_associated(plist%swigptr), .true.)

    ! Test a function that raises an exception
    TEST_THROW(call load_from_xml(plist, 'nonexistent_path.xml'))
    TEST_EQUALITY_CONST(c_associated(plist%swigptr), .true.)

    ! Get and set a vlaue
    call plist%set('myint', 4)
    call plist%get('myint', ival)
    TEST_EQUALITY(ival, 4)

    call plist%set('mydbl', 1.25_C_DOUBLE)
    call plist%get('mydbl', dval)
    TEST_FLOATING_EQUALITY(dval, 1.25_C_DOUBLE, epsilon(1.0_C_DOUBLE))

    bval = .false.
    call plist%set('mybool', true)
    call plist%get('mybool', bval)
    TEST_EQUALITY_CONST(bval, true)

    call plist%set('intarr', test_int)
    call plist%set('dblarr', test_dbl)

    TEST_EQUALITY(6, plist%get_length('intarr'))
    TEST_EQUALITY(4, plist%get_length('dblarr'))

    ! Wrong parameter type
    TEST_THROW(call plist%get('intarr', test_dbl))

    ! Wrong array size
    TEST_THROW(call plist%get('intarr', test_int(:4)))

    call plist%set('deleteme', 123)
    TEST_EQUALITY_CONST(plist%is_parameter('deleteme'), true)

    call plist%remove('deleteme')
    TEST_EQUALITY_CONST(plist%is_parameter('deleteme'), false)

    TEST_EQUALITY_CONST(c_associated(sublist%swigptr), .false.)
    sublist = plist%sublist('sublist')
    TEST_EQUALITY_CONST(c_associated(sublist%swigptr), .true.)

    call sublist%set('anotherval', 4.0d0)

    call sublist%set('stringval', 'some string!')
    TEST_EQUALITY(sublist%get_length('stringval'), 12)
    call sublist%get('stringval', sval)
    TEST_EQUALITY('some string!', sval)

    ! Set a string that's too long for sval
    TEST_NOTHROW(call sublist%set('stringval', 'the string is too damn long!'))
    TEST_THROW(call sublist%get('stringval', sval))

    call plist%set('sublist2', sublist)

    ! Add a sub-sublist
    sublist = sublist%sublist('subsublist')

    ! Add another parameter on the original sublist
    call sublist%set('late_arrival', 1.0d0)
    call sublist%release()

    ! XXX: make this interface cleaner. Currently getting a reference to a
    ! sublist requires that another parameterlist is allocated first
    call sublistref%create()
    call plist%get('sublist', sublistref)
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
