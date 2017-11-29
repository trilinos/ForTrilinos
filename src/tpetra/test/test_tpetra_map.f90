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

  type(TeuchosComm) :: comm
  integer(global_size_type), parameter :: invalid=-1
  integer(global_ordinal_type), parameter :: index_base=1

  SETUP_TEST()

#ifdef HAVE_MPI
  call comm%create(MPI_COMM_WORLD); CHECK_IERR()
#else
  call comm%create(); CHECK_IERR()
#endif

  ADD_SUBTEST_AND_RUN(TpetraMap_isOneToOne)
  ADD_SUBTEST_AND_RUN(TpetraMap_getGlobalNumElements)
  ADD_SUBTEST_AND_RUN(TpetraMap_getNodeNumElements)
  ADD_SUBTEST_AND_RUN(TpetraMap_getIndexBase)
  ADD_SUBTEST_AND_RUN(TpetraMap_getMinLocalIndex)
  ADD_SUBTEST_AND_RUN(TpetraMap_getMaxLocalIndex)
  ADD_SUBTEST_AND_RUN(TpetraMap_getMinGlobalIndex)
  ADD_SUBTEST_AND_RUN(TpetraMap_getMaxGlobalIndex)
  ADD_SUBTEST_AND_RUN(TpetraMap_getMinAllGlobalIndex)
  ADD_SUBTEST_AND_RUN(TpetraMap_getMaxAllGlobalIndex)
  ADD_SUBTEST_AND_RUN(TpetraMap_getLocalElement)
  ADD_SUBTEST_AND_RUN(TpetraMap_getGlobalElement)
  ADD_SUBTEST_AND_RUN(TpetraMap_getNodeElementList)
  ADD_SUBTEST_AND_RUN(TpetraMap_isNodeLocalElement)
  ADD_SUBTEST_AND_RUN(TpetraMap_isNodeGlobalElement)
  ADD_SUBTEST_AND_RUN(TpetraMap_isUniform)
  ADD_SUBTEST_AND_RUN(TpetraMap_isContiguous)
  ADD_SUBTEST_AND_RUN(TpetraMap_isDistributed)
  ADD_SUBTEST_AND_RUN(TpetraMap_isCompatible)
  ADD_SUBTEST_AND_RUN(TpetraMap_isSameAs)
  ADD_SUBTEST_AND_RUN(TpetraMap_locallySameAs)
  ADD_SUBTEST_AND_RUN(TpetraMap_getComm)
  ADD_SUBTEST_AND_RUN(TpetraMap_description)

  ! Methods listed below are marked as "may go away or change at any time"
  ADD_SUBTEST_AND_RUN(TpetraMap_removeEmptyProcesses)
  ADD_SUBTEST_AND_RUN(TpetraMap_replaceCommWithSubset)

  call comm%release()

  TEARDOWN_TEST()

contains

! ---------------------------------isOneToOne--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isOneToOne)
    logical(c_bool) :: bool
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global, indices(4)
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_ASSERT(Obj%isOneToOne())
    call Obj%release(); TEST_IERR()
    if (comm%getSize() > 1) then
      indices = [1, 2, 3, 4]
      call Obj%create(num_global, indices, index_base, comm); TEST_IERR()
      TEST_ASSERT((.not. Obj%isOneToOne()))
      call Obj%release(); TEST_IERR()
    end if
  END_FORTRILINOS_UNIT_TEST(TpetraMap_isOneToOne)

! ----------------------------getGlobalNumElements---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getGlobalNumElements)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_ASSERT(Obj%getGlobalNumElements() == num_global)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getGlobalNumElements)

! -----------------------------getNodeNumElements----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getNodeNumElements)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getNodeNumElements(), int(4, kind=size_type))
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getNodeNumElements)

! --------------------------------getIndexBase-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getIndexBase)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global, index_base_0
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getIndexBase(), int(1, kind=global_size_type))
    call Obj%release(); TEST_IERR()
    index_base_0 = 0
    call Obj%create(num_global, index_base_0, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getIndexBase(), int(0, kind=global_size_type))
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getIndexBase)

! ------------------------------getMinLocalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMinLocalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getMinLocalIndex(), int(1, kind=local_ordinal_type))
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMinLocalIndex)

! ------------------------------getMaxLocalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMaxLocalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getMaxLocalIndex(), int(4, kind=local_ordinal_type))
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMaxLocalIndex)

! -----------------------------getMinGlobalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMinGlobalIndex)
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: expected
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    expected = comm%getRank() * 4 + 1
    TEST_EQUALITY(Obj%getMinGlobalIndex(), expected)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMinGlobalIndex)

! -----------------------------getMaxGlobalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMaxGlobalIndex)
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: expected
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    expected = comm%getRank() * 4 + 4
    TEST_EQUALITY(Obj%getMaxGlobalIndex(), expected)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMaxGlobalIndex)

! ----------------------------getMinAllGlobalIndex---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMinAllGlobalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getMinAllGlobalIndex(), index_base)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMinAllGlobalIndex)

! ----------------------------getMaxAllGlobalIndex---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMaxAllGlobalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getMaxAllGlobalIndex(), num_global)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMaxAllGlobalIndex)

! ------------------------------getLocalElement------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getLocalElement)
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    globalindex = comm%getRank() * 4 + 1
    TEST_EQUALITY(Obj%getLocalElement(globalindex), 1)
    globalindex = comm%getRank() * 4 + 4
    TEST_EQUALITY(Obj%getLocalElement(globalindex), 4)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getLocalElement)

! ------------------------------getGlobalElement------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getGlobalElement)
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    integer(C_LONG_LONG) :: expected
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    localindex = 1
    expected = comm%getRank() * 4 + 1
    TEST_EQUALITY(Obj%getGlobalElement(localindex), expected)
    localindex = 4
    expected = comm%getRank() * 4 + 4
    TEST_EQUALITY(Obj%getGlobalElement(localindex), expected)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getGlobalElement)

! -----------------------------getNodeElementList----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getNodeElementList)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type), allocatable :: element_list(:)
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    allocate(element_list(4))
    ! TODO: The element list returned is junk
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    call Obj%getNodeElementList(element_list)
    deallocate(element_list)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getNodeElementList)

! -----------------------------isNodeLocalElement----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isNodeLocalElement)
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    localindex = 1
    TEST_ASSERT(Obj%isNodeLocalElement(localindex))
    localindex = 5
    TEST_ASSERT((.not. Obj%isNodeLocalElement(localindex)))
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_isNodeLocalElement)

! ----------------------------isNodeGlobalElement----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isNodeGlobalElement)
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()

    globalindex = int(comm%getRank() * 4 + 1, kind=global_ordinal_type)
    TEST_ASSERT(Obj%isNodeGlobalElement(globalindex))

    globalindex = num_global + 1
    TEST_ASSERT((.not. Obj%isNodeGlobalElement(globalindex)))

    call Obj%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isNodeGlobalElement)

! ---------------------------------isUniform---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isUniform)
    type(TpetraMap) :: Obj
    integer :: k
    integer(global_ordinal_type) :: num_global, elements(4)
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_ASSERT(Obj%isUniform())
    call Obj%release(); TEST_IERR()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm); TEST_IERR()
    TEST_ASSERT((.not. Obj%isUniform()))
    call Obj%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isUniform)

! --------------------------------isContiguous-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isContiguous)
    type(TpetraMap) :: Obj
    integer :: k
    integer(global_ordinal_type) :: num_global, elements(4)
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()
    TEST_ASSERT(Obj%isContiguous())
    call Obj%release(); TEST_IERR()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm); TEST_IERR()
    TEST_ASSERT((.not. Obj%isContiguous()))
    call Obj%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isContiguous)

! -------------------------------isDistributed-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isDistributed)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global, elements(4)
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()
    if (comm%getSize() == 1) then
      TEST_ASSERT((.not. Obj%isDistributed()))
    else
      TEST_ASSERT(Obj%isDistributed())
    end if
    call Obj%release(); TEST_IERR()

    ! All elements have 4 entries and map has only 4 entries so the map should
    ! not be distributed
    num_global = 4
    elements = [1, 2, 3, 4]
    call Obj%create(num_global, elements, index_base, comm); TEST_IERR()
    TEST_ASSERT((.not. Obj%isDistributed()))
    call Obj%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isDistributed)

! --------------------------------isCompatible-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isCompatible)
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    integer(global_ordinal_type) :: num_global, elements(4)
    integer :: k
    num_global = 4*comm%getSize()

    call Obj1%create(num_global, index_base, comm); TEST_IERR()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj2%create(num_global, elements, index_base, comm); TEST_IERR()

    ! Cyclic map should be compatible
    TEST_ASSERT(Obj1%isCompatible(Obj2))

    call Obj2%release(); TEST_IERR()

    num_global = num_global + 10
    call Obj2%create(num_global, index_base, comm); TEST_IERR()

    TEST_ASSERT((.not. Obj1%isCompatible(Obj2)))

    call Obj2%release(); TEST_IERR()
    call Obj1%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isCompatible)

! ----------------------------------isSameAs---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isSameAs)
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, elements(4)
    integer :: k
    num_global = 4*comm%getSize()

    call Obj1%create(num_global, index_base, comm); TEST_IERR()

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      call Obj2%create(num_global, elements, index_base, comm); TEST_IERR()

      ! Cyclic map should not be SameAs
      TEST_ASSERT((.not. Obj1%isSameAs(Obj2)))

      call Obj2%release(); TEST_IERR()

    end if

    num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm); TEST_IERR()

    TEST_ASSERT(Obj1%isSameAs(Obj2))

    call Obj1%release(); TEST_IERR()
    call Obj2%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isSameAs)

! -------------------------------locallySameAs-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_locallySameAs)
    type(TpetraMap) :: Obj1, Obj2
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, elements(4)
    integer :: k
    num_global = 4*comm%getSize()

    call Obj1%create(num_global, index_base, comm); TEST_IERR()

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      call Obj2%create(num_global, elements, index_base, comm); TEST_IERR()

      ! Cyclic map should not be locallySameAs
      TEST_ASSERT((.not. Obj1%locallySameAs(Obj2)))

      call Obj2%release(); TEST_IERR()

    end if

    num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm); TEST_IERR()
    TEST_ASSERT(Obj1%locallySameAs(Obj2))

    call Obj1%release(); TEST_IERR()
    call Obj2%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_locallySameAs)

! ----------------------------------getComm----------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getComm)
    type(TpetraMap) :: Obj
    type(TeuchosComm) :: tcomm
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    ! We only test the comm returned has the same rank and size.  More
    ! comprehensive testing is (assumed to be) done in Teuchos itself.
    tcomm = Obj%getComm(); TEST_IERR()
    TEST_ASSERT(tcomm%getRank() == comm%getRank())
    TEST_ASSERT(tcomm%getSize() == comm%getSize())

    call Obj%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getComm)

! --------------------------------description--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_description)
    type(TpetraMap) :: Obj
    type(string) :: desciption
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    desciption = Obj%description(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_description)

! ----------------------------removeEmptyProcesses---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_removeEmptyProcesses)
    ! TODO: Implement this test?
    if (comm%getRank()==0) &
      write(*,*) 'removeEmptyProcesses: Test is not yet Implemented'
  END_FORTRILINOS_UNIT_TEST(TpetraMap_removeEmptyProcesses)

! ---------------------------replaceCommWithSubset---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_replaceCommWithSubset)
    ! TODO: Implement this test?
    if (comm%getRank()==0) &
      write(*,*) 'replaceCommWithSubset: Test is not yet implemented'
  END_FORTRILINOS_UNIT_TEST(TpetraMap_replaceCommWithSubset)

end program test_TpetraMap
