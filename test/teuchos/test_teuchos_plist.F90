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
  character(len=26), parameter :: FILENAME='test_teuchos_plist.F90'
  integer, parameter :: dp = kind(0.d0)

  SETUP_TEST()

  ADD_SUBTEST_AND_RUN(TeuchosPList_Basic)

  ADD_SUBTEST_AND_RUN(ParameterList_print)
  ADD_SUBTEST_AND_RUN(ParameterList_remove)
  ADD_SUBTEST_AND_RUN(ParameterList_is_parameter)
  ADD_SUBTEST_AND_RUN(ParameterList_sublist)
  ADD_SUBTEST_AND_RUN(ParameterList_get_real)
  ADD_SUBTEST_AND_RUN(ParameterList_get_integer)
  ADD_SUBTEST_AND_RUN(ParameterList_get_longlong)
  ADD_SUBTEST_AND_RUN(ParameterList_get_logical)
  ADD_SUBTEST_AND_RUN(ParameterList_get_string)
  ADD_SUBTEST_AND_RUN(ParameterList_get_arr_real)
  ADD_SUBTEST_AND_RUN(ParameterList_get_arr_integer)
  ADD_SUBTEST_AND_RUN(ParameterList_get_arr_longlong)

  TEARDOWN_TEST()
contains

  FORTRILINOS_UNIT_TEST(TeuchosPList_Basic)

    type(ParameterList) :: plist, sublist, sublistref
    integer(C_INT), dimension(6) :: test_int = (/ -1, 1, 3, 3, 5, 7 /)
    real(C_DOUBLE), dimension(4) :: test_dbl = (/ 0.1d0, 1.9d0, -2.0d0, 4.0d0 /)

    integer(C_INT), pointer :: iarr(:) => NULL()
    real(C_DOUBLE), pointer :: darr(:) => NULL()
    integer(C_INT) :: ival
    real(C_DOUBLE) :: dval
    logical :: bval
    character(kind=C_CHAR, len=:), allocatable :: sval

    OUT0('Starting TeuchosPList_Basic!')

    plist = ParameterList('myname'); TEST_IERR()
    TEST_ASSERT(c_associated(plist%swigdata%cptr))

    ! Test a function that raises an exception
    TEST_THROW(call load_from_xml(plist, 'nonexistent_path.xml'))
    TEST_ASSERT(c_associated(plist%swigdata%cptr))

    ! Get and set a value
    call plist%set('myint', 4)
    ival = plist%get_integer('myint')
    TEST_EQUALITY(ival, 4)

    call plist%set('mydbl', 1.25_C_DOUBLE)
    dval = plist%get_real('mydbl')
    TEST_FLOATING_EQUALITY(dval, 1.25_C_DOUBLE, epsilon(1.0_C_DOUBLE))

    bval = .false.
    call plist%set('mybool', .true.)
    bval = plist%get_logical('mybool')
    TEST_ASSERT(bval)

    call plist%set('intarr', test_int)
    call plist%set('dblarr', test_dbl)

    ! FIXME: without these allocation, the test segfaults when compiled with flang
    ! According to [this](https://stackoverflow.com/questions/42140832/automatic-array-allocation-upon-assignment-in-fortran) it
    ! seems that it should be supported in F2003. So it is not clear why flang fails.
    iarr => plist%get_arr_integer('intarr')
    darr => plist%get_arr_real('dblarr')

    ! Wrong parameter type
    TEST_THROW(iarr => plist%get_arr_integer('dblarr'))

    call plist%set('deleteme', 123)
    TEST_ASSERT(plist%is_parameter('deleteme'))

    call plist%remove('deleteme')
    TEST_ASSERT((.not. plist%is_parameter('deleteme')))

    TEST_ASSERT((.not. c_associated(sublist%swigdata%cptr)))
    sublist = plist%sublist('sublist')
    TEST_ASSERT(c_associated(sublist%swigdata%cptr))

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

  ! ----------------------------------print----------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_print)
    type(ParameterList) :: list
    OUT0("Starting ParameterList_print!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('someint', 1)
    call list%set('somedouble', 2.5_dp)
    call list%set('somestring', 'Hello World')
    call list%print(); TEST_IERR()

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_print!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_print)

  ! ----------------------------------remove---------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_remove)
    type(ParameterList) :: list
    OUT0("Starting ParameterList_remove!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('test', 0)
    TEST_ASSERT(list%is_parameter('test'))

    call list%remove('test')
    TEST_ASSERT(.not. list%is_parameter('test'))

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_remove!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_remove)

  ! -------------------------------is_parameter------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_is_parameter)
    type(ParameterList) :: list
    OUT0("Starting ParameterList_is_parameter!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('p1', 1)
    call list%set('p2', 2)
    TEST_ASSERT(list%is_parameter('p1') .and. list%is_parameter('p2'))
    TEST_ASSERT(.not. list%is_parameter('p3'))

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_is_parameter!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_is_parameter)

  ! ---------------------------------sublist---------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_sublist)
    type(ParameterList) :: list, sub
    OUT0("Starting ParameterList_sublist!")

    list = ParameterList('mylist'); TEST_IERR()
    sub = list%sublist('mysub'); TEST_IERR()
    call sub%set('myint', 0)
    TEST_ASSERT(sub%is_parameter('myint'))

    call sub%release(); TEST_IERR()
    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_sublist!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_sublist)

  ! ---------------------------------get_real--------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_get_real)
    type(ParameterList) :: list
    real(dp) :: var
    OUT0("Starting ParameterList_get_real!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('myvar', 1.0_dp); TEST_IERR()
    var = list%get_real('myvar')
    TEST_FLOATING_EQUALITY(1.0_dp, var, epsilon(1.0_dp))

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_get_real!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_get_real)

  ! -------------------------------get_integer-------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_get_integer)
    type(ParameterList) :: list
    integer :: var
    OUT0("Starting ParameterList_get_integer!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('myvar', 1); TEST_IERR()
    var = list%get_integer('myvar')
    TEST_EQUALITY(1, var)

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_get_integer!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_get_integer)

  ! -------------------------------get_longlong------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_get_longlong)
    type(ParameterList) :: list
    integer, parameter :: intk = selected_int_kind(18)
    integer(intk) :: var
    OUT0("Starting ParameterList_get_longlong!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('myvar', 1_intk); TEST_IERR()
    var = list%get_longlong('myvar')
    TEST_EQUALITY(1_intk, var)

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_get_longlong!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_get_longlong)

  ! -------------------------------get_logical-------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_get_logical)
    type(ParameterList) :: list
    logical :: var
    OUT0("Starting ParameterList_get_logical!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('myvar', .true.); TEST_IERR()
    var = list%get_logical('myvar')
    TEST_ASSERT(var)

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_get_logical!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_get_logical)

  ! --------------------------------get_string-------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_get_string)
    type(ParameterList) :: list
    character(kind=C_CHAR, len=:), allocatable :: var
    OUT0("Starting ParameterList_get_string!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('myvar', 'Test'); TEST_IERR()
    var = list%get_string('myvar')
    TEST_EQUALITY(var, 'Test')

    deallocate(var)
    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_get_string!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_get_string)

  ! -------------------------------get_arr_real------------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_get_arr_real)
    type(ParameterList) :: list
    real(dp), pointer :: var(:) => NULL()
    real(dp), dimension(3) :: myvec = (/ -3.2_dp, 0.2_dp, 5.8_dp /)
    OUT0("Starting ParameterList_get_arr_real!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('rarr', myvec); TEST_IERR()
    var => list%get_arr_real('rarr')
    TEST_FLOATING_ARRAY_EQUALITY(var, myvec, epsilon(var(1)))

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_get_arr_real!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_get_arr_real)

  ! -----------------------------get_arr_integer------------------------------ !
  FORTRILINOS_UNIT_TEST(ParameterList_get_arr_integer)
    type(ParameterList) :: list
    integer, pointer :: var(:) => NULL()
    integer, dimension(3) :: myvec = (/ -3, 0, 6 /)
    OUT0("Starting ParameterList_get_arr_integer!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('iarr', myvec); TEST_IERR()
    var => list%get_arr_integer('iarr')
    TEST_ARRAY_EQUALITY(var, myvec)

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_get_arr_integer!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_get_arr_integer)

  ! -----------------------------get_arr_longlong----------------------------- !
  FORTRILINOS_UNIT_TEST(ParameterList_get_arr_longlong)
    type(ParameterList) :: list
    integer, parameter :: intk = selected_int_kind(18)
    integer(intk), pointer :: var(:) => NULL()
    integer(intk), dimension(3) :: myvec = (/ -3, 0, 6 /)
    OUT0("Starting ParameterList_get_arr_longlong!")

    list = ParameterList('mylist'); TEST_IERR()
    call list%set('larr', myvec); TEST_IERR()
    var => list%get_arr_longlong('larr')
    TEST_ARRAY_EQUALITY(var, myvec)

    call list%release(); TEST_IERR()

    OUT0("Finished ParameterList_get_arr_longlong!")

  END_FORTRILINOS_UNIT_TEST(ParameterList_get_arr_longlong)

end program test_TeuchosPList
