! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program test_TpetraMap
#include "ForTrilinosTpetra_config.hpp"
#include "FortranTestMacros.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  DECLARE_TEST_VARIABLES()
  type(TeuchosComm) :: comm

  INITIALIZE_TEST()

#ifdef HAVE_MPI
  call comm%create(MPI_COMM_WORLD)
  CHECK_IERR()
#else
  call comm%create()
  CHECK_IERR()
#endif

  ADD_SUBTEST_AND_RUN(test_isOneToOne)
  ADD_SUBTEST_AND_RUN(test_getGlobalNumElements)
  ADD_SUBTEST_AND_RUN(test_getNodeNumElements)
  ADD_SUBTEST_AND_RUN(test_getIndexBase)
  ADD_SUBTEST_AND_RUN(test_getMinLocalIndex)
  ADD_SUBTEST_AND_RUN(test_getMaxLocalIndex)
  ADD_SUBTEST_AND_RUN(test_getMinGlobalIndex)
  ADD_SUBTEST_AND_RUN(test_getMaxGlobalIndex)
  ADD_SUBTEST_AND_RUN(test_getMinAllGlobalIndex)
  ADD_SUBTEST_AND_RUN(test_getMaxAllGlobalIndex)
  ADD_SUBTEST_AND_RUN(test_getLocalElement)
  ADD_SUBTEST_AND_RUN(test_getGlobalElement)
  ADD_SUBTEST_AND_RUN(test_getNodeElementList)
  ADD_SUBTEST_AND_RUN(test_isNodeLocalElement)
  ADD_SUBTEST_AND_RUN(test_isNodeGlobalElement)
  ADD_SUBTEST_AND_RUN(test_isUniform)
  ADD_SUBTEST_AND_RUN(test_isContiguous)
  ADD_SUBTEST_AND_RUN(test_isDistributed)
  ADD_SUBTEST_AND_RUN(test_isCompatible)
  ADD_SUBTEST_AND_RUN(test_isSameAs)
  ADD_SUBTEST_AND_RUN(test_locallySameAs)
  ADD_SUBTEST_AND_RUN(test_getComm)
  ADD_SUBTEST_AND_RUN(test_description)
  ADD_SUBTEST_AND_RUN(test_removeEmptyProcesses)
  ADD_SUBTEST_AND_RUN(test_replaceCommWithSubset)

  call comm%release()
  CHECK_IERR()

  SHUTDOWN_TEST()

contains

! ---------------------------------isOneToOne--------------------------------- !
  integer function test_isOneToOne()
    integer :: jerr
    logical(c_bool) :: bool
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global, index_base, indices(4)
    jerr = 0
    index_base = 1
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_isOneToOne)

    bool = Obj%isOneToOne(); TEST_FOR_IERR(test_isOneToOne)
    if (.not. bool) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isOneToOne: Expected map to be one to one"
    end if

    call Obj%release(); TEST_FOR_IERR(test_isOneToOne)

    if (comm%getSize() > 1) then
      indices = [1, 2, 3, 4]
      call Obj%create(num_global, indices, index_base, comm)
      TEST_FOR_IERR(test_isOneToOne)

      bool = Obj%isOneToOne(); TEST_FOR_IERR(test_isOneToOne)
      if (bool) then
        jerr = jerr + 1
        if (comm%getRank() == 0) &
          write(*,*) "isOneToOne: Expected map to NOT be one to one"
      end if

      call Obj%release(); TEST_FOR_IERR(test_isOneToOne)

    end if

    SET_ERROR_COUNT_AND_RETURN(test_isOneToOne, jerr)

  end function

! ----------------------------getGlobalNumElements---------------------------- !
  integer function test_getGlobalNumElements()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    index_base = 1
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getGlobalNumElements)

    fresult = Obj%getGlobalNumElements()
    TEST_FOR_IERR(test_getGlobalNumElements)
    if (fresult /= num_global) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getGlobalNumElements: Expected ", num_global, " elements, got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getGlobalNumElements)

    SET_ERROR_COUNT_AND_RETURN(test_getGlobalNumElements, jerr)

  end function

! -----------------------------getNodeNumElements----------------------------- !
  integer function test_getNodeNumElements()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_SIZE_T) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    index_base = 1
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getNodeNumElements)

    fresult = Obj%getNodeNumElements()
    if (fresult /= 4) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getNodeNumElements: Expected ", 4, " elements, got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getNodeNumElements)

    SET_ERROR_COUNT_AND_RETURN(test_getNodeNumElements, jerr)

  end function

! --------------------------------getIndexBase-------------------------------- !
  integer function test_getIndexBase()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getIndexBase)

    fresult = Obj%getIndexBase()
    TEST_FOR_IERR(test_getIndexBase)
    if (fresult /= 1) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getIndexBase: Expected indexBase = ", 1, " got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getIndexBase)

    index_base = 0
    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getIndexBase)

    fresult = Obj%getIndexBase()
    TEST_FOR_IERR(test_getIndexBase)
    if (fresult /= 0) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getIndexBase: Expected indexBase = ", 0, " got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getIndexBase)

    SET_ERROR_COUNT_AND_RETURN(test_getIndexBase, jerr)

  end function

! ------------------------------getMinLocalIndex------------------------------ !
  integer function test_getMinLocalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getMinLocalIndex)

    fresult = Obj%getMinLocalIndex()
    TEST_FOR_IERR(test_getMinLocalIndex)
    if (fresult /= 1) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMinLocalIndex: Expected minLocalIndex = ", 1, " got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getMinLocalIndex)

    SET_ERROR_COUNT_AND_RETURN(test_getMinLocalIndex, jerr)

  end function

! ------------------------------getMaxLocalIndex------------------------------ !
  integer function test_getMaxLocalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getMaxLocalIndex)

    fresult = Obj%getMaxLocalIndex()
    TEST_FOR_IERR(test_getMaxLocalIndex)
    if (fresult /= 4) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMaxLocalIndex: Expected maxLocalIndex = ", 4, " got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getMaxLocalIndex)

    SET_ERROR_COUNT_AND_RETURN(test_getMaxLocalIndex, jerr)

  end function

! -----------------------------getMinGlobalIndex------------------------------ !
  integer function test_getMinGlobalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getMinGlobalIndex)

    fresult = Obj%getMinGlobalIndex()
    TEST_FOR_IERR(test_getMinGlobalIndex)
    expected = comm%getRank() * 4 + 1
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMindGlobalIndex: Expected minGlobalIndex = ", expected, " got ", fresult
    end if
    call Obj%release()
    TEST_FOR_IERR(test_getMinGlobalIndex)

    SET_ERROR_COUNT_AND_RETURN(test_getMinGlobalIndex, jerr)

  end function

! -----------------------------getMaxGlobalIndex------------------------------ !
  integer function test_getMaxGlobalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getMaxGlobalIndex)

    fresult = Obj%getMaxGlobalIndex()
    TEST_FOR_IERR(test_getMaxGlobalIndex)
    expected = comm%getRank() * 4 + 4
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMaxGlobalIndex: Expected maxGloblIndex = ", expected, " got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getMaxGlobalIndex)

    SET_ERROR_COUNT_AND_RETURN(test_getMaxGlobalIndex, jerr)

  end function

! ----------------------------getMinAllGlobalIndex---------------------------- !
  integer function test_getMinAllGlobalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getMinAllGlobalIndex)

    fresult = Obj%getMinAllGlobalIndex()
    TEST_FOR_IERR(test_getMinAllGlobalIndex)
    if (fresult /= index_base) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMinAllGlobalIndex: Expected minAllGlobalIndex = ", index_base, " got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getMinAllGlobalIndex)

    SET_ERROR_COUNT_AND_RETURN(test_getMinAllGlobalIndex, jerr)

  end function

! ----------------------------getMaxAllGlobalIndex---------------------------- !
  integer function test_getMaxAllGlobalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getMaxAllGlobalIndex)

    fresult = Obj%getMaxAllGlobalIndex()
    TEST_FOR_IERR(test_getMaxAllGlobalIndex)
    if (fresult /= num_global) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMaxAllGlobalIndex: Expected maxAllGlobalIndex = ", num_global, " got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getMaxAllGlobalIndex)

    SET_ERROR_COUNT_AND_RETURN(test_getMaxAllGlobalIndex, jerr)

  end function

! ------------------------------getLocalElement------------------------------- !
  integer function test_getLocalElement()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getLocalElement)

    globalindex = comm%getRank() * 4 + 1
    fresult = Obj%getLocalElement(globalindex)
    TEST_FOR_IERR(test_getLocalElement)
    if (fresult /= 1) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getLocalElement: Expected local element = ", 1, " got ", fresult
    end if

    globalindex = comm%getRank() * 4 + 4
    fresult = Obj%getLocalElement(globalindex)
    TEST_FOR_IERR(test_getLocalElement)
    if (fresult /= 4) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getLocalElement: Expected local element = ", 4, " got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getLocalElement)

    SET_ERROR_COUNT_AND_RETURN(test_getLocalElement, jerr)

  end function

! ------------------------------getGlobalElement------------------------------ !
  integer function test_getGlobalElement()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getGlobalElement)

    localindex = 1
    expected = comm%getRank() * 4 + 1
    fresult = Obj%getGlobalElement(localindex)
    TEST_FOR_IERR(test_getGlobalElement)
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getGlobalElement: Expected global element = ", expected, " got ", fresult
    end if

    localindex = 4
    expected = comm%getRank() * 4 + 4
    fresult = Obj%getGlobalElement(localindex)
    TEST_FOR_IERR(test_getGlobalElement)
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getGlobalElement: Expected global element = ", expected, " got ", fresult
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getGlobalElement)

    SET_ERROR_COUNT_AND_RETURN(test_getGlobalElement, jerr)

  end function

! -----------------------------getNodeElementList----------------------------- !
  integer function test_getNodeElementList()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(local_ordinal_type), parameter :: num_local = 4
    integer(global_ordinal_type) :: fresult(num_local)
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = num_local*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getNodeElementList)

    call Obj%getNodeElementList(fresult)
    TEST_FOR_IERR(test_getNodeElementList)

    ! TODO: getNodeElementList should be modified to be a subroutine that takes
    ! TODO: the element list as a return argument.  Otherwise, there is no
    ! TODO: (current) way to get to the data in the ArrayView returned.
    if (comm%getRank() == 0) &
      write(*,*) "getNodeElementList: Test not yet implemented"

    TEST_FOR_IERR(test_getNodeElementList)

    call Obj%release()
    TEST_FOR_IERR(test_getNodeElementList)

    SET_ERROR_COUNT_AND_RETURN(test_getNodeElementList, jerr)

  end function

! -----------------------------isNodeLocalElement----------------------------- !
  integer function test_isNodeLocalElement()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_isNodeLocalElement)

    localindex = 1
    fresult = Obj%isNodeLocalElement(localindex)
    TEST_FOR_IERR(test_isNodeLocalElement)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeLocalElement: expected 1 to be a local index"
    end if

    localindex = 5
    fresult = Obj%isNodeLocalElement(localindex)
    TEST_FOR_IERR(test_isNodeLocalElement)
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeLocalElement: expected 5 to NOT be a local index"
    end if

    call Obj%release()
    TEST_FOR_IERR(test_isNodeLocalElement)

    SET_ERROR_COUNT_AND_RETURN(test_isNodeLocalElement, jerr)

  end function

! ----------------------------isNodeGlobalElement----------------------------- !
  integer function test_isNodeGlobalElement()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_isNodeGlobalElement)

    globalindex = int(comm%getRank() * 4 + 1, kind=global_ordinal_type)
    fresult = Obj%isNodeGlobalElement(globalindex)
    TEST_FOR_IERR(test_isNodeGlobalElement)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeGlobalElement: expected 1 to be a global index"
    end if
    globalindex = num_global + 1

    fresult = Obj%isNodeGlobalElement(globalindex)
    TEST_FOR_IERR(test_isNodeGlobalElement)
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeGlobalElement: expected 5 to NOT be a global index"
    end if

    call Obj%release()
    TEST_FOR_IERR(test_isNodeGlobalElement)

    SET_ERROR_COUNT_AND_RETURN(test_isNodeGlobalElement, jerr)

  end function

! ---------------------------------isUniform---------------------------------- !
  integer function test_isUniform()
    integer :: jerr
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer :: k
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_isUniform)

    fresult = Obj%isUniform()
    TEST_FOR_IERR(test_isUniform)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isUniform: expected map to be uniform"
    end if

    call Obj%release()
    TEST_FOR_IERR(test_isUniform)

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm);
    TEST_FOR_IERR(test_isUniform)

    fresult = Obj%isUniform()
    TEST_FOR_IERR(test_isUniform)
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isUniform: expected map to NOT be uniform"
    end if

    call Obj%release()
    TEST_FOR_IERR(test_isUniform)

    SET_ERROR_COUNT_AND_RETURN(test_isUniform, jerr)

  end function

! --------------------------------isContiguous-------------------------------- !
  integer function test_isContiguous()
    integer :: jerr
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer :: k
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_isContiguous)

    fresult = Obj%isContiguous()
    TEST_FOR_IERR(test_isContiguous)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isContiguous: expected map to be uniform"
    end if

    call Obj%release()
    TEST_FOR_IERR(test_isContiguous)

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm);
    TEST_FOR_IERR(test_isContiguous)

    fresult = Obj%isContiguous()
    TEST_FOR_IERR(test_isContiguous)
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isContiguous: expected map to NOT be uniform"
    end if

    call Obj%release()
    TEST_FOR_IERR(test_isContiguous)

    SET_ERROR_COUNT_AND_RETURN(test_isContiguous, jerr)

  end function

! -------------------------------isDistributed-------------------------------- !
  integer function test_isDistributed()
    integer :: jerr
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_isDistributed)

    fresult = Obj%isDistributed()
    TEST_FOR_IERR(test_isDistributed)
    if (comm%getSize() == 1) then

      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isDistributed: expected map to NOT be distributed"
      end if

    else

      if (.not. fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isDistributed: expected map to be distributed"
      end if

    end if

    call Obj%release()
    TEST_FOR_IERR(test_isDistributed)

    ! All elements have 4 entries and map has only 4 entries so the map should
    ! not be distributed
    num_global = 4
    elements = [1, 2, 3, 4]
    call Obj%create(num_global, elements, index_base, comm)
    TEST_FOR_IERR(test_isDistributed)

    fresult = Obj%isDistributed()
    TEST_FOR_IERR(test_isDistributed)
    if (comm%getSize() == 1) then

      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isDistributed: expected map to NOT be distributed"
      end if

    else

      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isDistributed: expected map to NOT be distributed"
      end if

    end if

    call Obj%release()
    TEST_FOR_IERR(test_isDistributed)

    SET_ERROR_COUNT_AND_RETURN(test_isDistributed, jerr)

  end function

! --------------------------------isCompatible-------------------------------- !
  integer function test_isCompatible()
    integer :: jerr
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    integer :: k
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj1%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_isCompatible)

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj2%create(num_global, elements, index_base, comm);
    TEST_FOR_IERR(test_isCompatible)

    ! Cyclic map should be compatible
    fresult = Obj1%isCompatible(Obj2)
    TEST_FOR_IERR(test_isCompatible)
    if (.not. fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isCompatible: Expected maps to be compatible"
    end if

    call Obj2%release()
    TEST_FOR_IERR(test_isCompatible)

    num_global = num_global + 10
    call Obj2%create(num_global, index_base, comm);
    TEST_FOR_IERR(test_isCompatible)

    fresult = Obj1%isCompatible(Obj2)
    TEST_FOR_IERR(test_isCompatible)
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "isCompatible: Expected maps to NOT be compatible"
    end if

    call Obj2%release()
    TEST_FOR_IERR(test_isCompatible)

    call Obj1%release()
    TEST_FOR_IERR(test_isCompatible)

    SET_ERROR_COUNT_AND_RETURN(test_isCompatible, jerr)

  end function

! ----------------------------------isSameAs---------------------------------- !
  integer function test_isSameAs()
    integer :: jerr
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    logical(C_BOOL) :: fresult
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    integer(global_size_type) :: invalid
    integer :: k
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj1%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_isSameAs)

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      call Obj2%create(num_global, elements, index_base, comm);
      TEST_FOR_IERR(test_isSameAs)

      ! Cyclic map should not be SameAs
      fresult = Obj1%isSameAs(Obj2)
      TEST_FOR_IERR(test_isSameAs)
      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isSameAs: Expected maps to NOT be same"
      end if

      call Obj2%release()
      TEST_FOR_IERR(test_isSameAs)

    end if

    invalid = -1; num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm)
    TEST_FOR_IERR(test_isSameAs)

    fresult = Obj1%isSameAs(Obj2)
    TEST_FOR_IERR(test_isSameAs)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "isSameAs: Expected maps to be same"
    end if

    call Obj1%release()
    TEST_FOR_IERR(test_isSameAs)

    call Obj2%release()
    TEST_FOR_IERR(test_isSameAs)

    SET_ERROR_COUNT_AND_RETURN(test_isSameAs, jerr)

  end function

! -------------------------------locallySameAs-------------------------------- !
  integer function test_locallySameAs()
    integer :: jerr
    type(TpetraMap) :: Obj1, Obj2
    logical(C_BOOL) :: fresult
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    integer(global_size_type) :: invalid
    integer :: k
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj1%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_locallySameAs)

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      call Obj2%create(num_global, elements, index_base, comm);
      TEST_FOR_IERR(test_locallySameAs)

      ! Cyclic map should not be locallySameAs
      fresult = Obj1%locallySameAs(Obj2)
      TEST_FOR_IERR(test_locallySameAs)
      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "locallySameAs: Expected maps to NOT be same"
      end if

      call Obj2%release()
      TEST_FOR_IERR(test_locallySameAs)

    end if

    invalid = -1; num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm)
    TEST_FOR_IERR(test_locallySameAs)

    fresult = Obj1%locallySameAs(Obj2)
    TEST_FOR_IERR(test_locallySameAs)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "locallySameAs: Expected maps to be same"
    end if

    call Obj1%release()
    TEST_FOR_IERR(test_locallySameAs)

    call Obj2%release()
    TEST_FOR_IERR(test_locallySameAs)

    SET_ERROR_COUNT_AND_RETURN(test_locallySameAs, jerr)

  end function

! ----------------------------------getComm----------------------------------- !
  integer function test_getComm()
    integer :: jerr
    type(TpetraMap) :: Obj
    type(TeuchosComm) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_getComm)

    ! We only test the comm returned has the same rank and size.  More
    ! comprehensive testing is (assumed to be) done in Teuchos itself.
    fresult = Obj%getComm()
    TEST_FOR_IERR(test_getComm)
    if (fresult%getRank() /= comm%getRank()) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "getComm: expected ranks to be same"
    end if

    if (fresult%getSize() /= comm%getSize()) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "getComm: expected sizes to be same"
    end if

    call Obj%release()
    TEST_FOR_IERR(test_getComm)

    SET_ERROR_COUNT_AND_RETURN(test_getComm, jerr)

  end function

! --------------------------------description--------------------------------- !
  integer function test_description()
    integer :: jerr
    type(TpetraMap) :: Obj
    type(string) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    TEST_FOR_IERR(test_description)
    fresult = Obj%description()
    TEST_FOR_IERR(test_description)
    SET_ERROR_COUNT_AND_RETURN(test_description, jerr)
  end function

! ----------------------------removeEmptyProcesses---------------------------- !
  integer function test_removeEmptyProcesses()
    integer :: jerr
    jerr = 0
    ! TODO: Implement this test?
    if (comm%getRank()==0) &
      write(*,*) 'removeEmptyProcesses: Test is not yet Implemented'
    SET_ERROR_COUNT_AND_RETURN(test_removeEmptyProcesses, jerr)
  end function

! ---------------------------replaceCommWithSubset---------------------------- !
  integer function test_replaceCommWithSubset()
    integer :: jerr
    jerr = 0
    ! TODO: Implement this test?
    if (comm%getRank()==0) &
      write(*,*) 'replaceCommWithSubset: Test is not yet implemented'
    SET_ERROR_COUNT_AND_RETURN(test_replaceCommWithSubset, jerr)
  end function

end program test_TpetraMap
