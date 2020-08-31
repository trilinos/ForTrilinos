! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program test_TpetraMap
#include "ForTrilinos_config.h"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  implicit none
  type(TeuchosComm) :: comm
  character(len=256), parameter :: FILENAME="test_tpetra_map.F90"

  SETUP_TEST()

#if FORTRILINOS_USE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm(); FORTRILINOS_CHECK_IERR()
#endif

  ADD_SUBTEST_AND_RUN(TpetraMap_isOneToOne)
  ADD_SUBTEST_AND_RUN(TpetraMap_getGlobalNumElements)
  ADD_SUBTEST_AND_RUN(TpetraMap_getNodeNumElements)
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
  ADD_SUBTEST_AND_RUN(TpetraMap_isLocallyFitted)

  ! Methods listed below are marked as "may go away or change at any time"
  !ADD_SUBTEST_AND_RUN(TpetraMap_removeEmptyProcesses)
  !ADD_SUBTEST_AND_RUN(TpetraMap_replaceCommWithSubset)

  call comm%release();  TEST_IERR()

  TEARDOWN_TEST()

contains

! ---------------------------------isOneToOne--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isOneToOne)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global, indices(4)
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    TEST_ASSERT(Obj%isOneToOne())
    call Obj%release(); TEST_IERR()
    if (comm%getSize() > 1) then
      indices = [1, 2, 3, 4]
      Obj = TpetraMap(num_global, indices, comm); TEST_IERR()
      TEST_ASSERT((.not. Obj%isOneToOne()))
      call Obj%release(); TEST_IERR()
    end if
  END_FORTRILINOS_UNIT_TEST(TpetraMap_isOneToOne)

! ----------------------------getGlobalNumElements---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getGlobalNumElements)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    TEST_ASSERT(Obj%getGlobalNumElements() == num_global)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getGlobalNumElements)

! -----------------------------getNodeNumElements----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getNodeNumElements)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getNodeNumElements(), 4)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getNodeNumElements)

! ------------------------------getMinLocalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMinLocalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getMinLocalIndex(), 1)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMinLocalIndex)

! ------------------------------getMaxLocalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMaxLocalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getMaxLocalIndex(), 4)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMaxLocalIndex)

! -----------------------------getMinGlobalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMinGlobalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: expected
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    expected = comm%getRank() * 4 + 1
    TEST_EQUALITY(Obj%getMinGlobalIndex(), expected)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMinGlobalIndex)

! -----------------------------getMaxGlobalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMaxGlobalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: expected
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    expected = comm%getRank() * 4 + 4
    TEST_EQUALITY(Obj%getMaxGlobalIndex(), expected)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMaxGlobalIndex)

! ----------------------------getMinAllGlobalIndex---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMinAllGlobalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global, one
    num_global = 4*comm%getSize()
    one = 1
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getMinAllGlobalIndex(), one)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMinAllGlobalIndex)

! ----------------------------getMaxAllGlobalIndex---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMaxAllGlobalIndex)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    TEST_EQUALITY(Obj%getMaxAllGlobalIndex(), num_global)
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMaxAllGlobalIndex)

! ------------------------------getLocalElement------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getLocalElement)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: globalindex
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
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
    integer(global_ordinal_type) :: expected
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
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
    integer(global_ordinal_type), pointer :: element_list(:)
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    element_list => Obj%getNodeElementList()
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_getNodeElementList)

! -----------------------------isNodeLocalElement----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isNodeLocalElement)
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    localindex = 1
    TEST_ASSERT(Obj%isNodeLocalElement(localindex))
    localindex = 5
    TEST_ASSERT((.not. Obj%isNodeLocalElement(localindex)))
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_isNodeLocalElement)

! ----------------------------isNodeGlobalElement----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isNodeGlobalElement)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: globalindex
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()

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

    Obj = TpetraMap(num_global, comm); TEST_IERR()
    TEST_ASSERT(Obj%isUniform())
    call Obj%release(); TEST_IERR()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    Obj = TpetraMap(num_global, elements, comm); TEST_IERR()
    TEST_ASSERT((.not. Obj%isUniform()))
    call Obj%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isUniform)

! --------------------------------isContiguous-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isContiguous)
    type(TpetraMap) :: Obj
    integer :: k
    integer(global_ordinal_type) :: num_global, elements(4)
    num_global = 4*comm%getSize()

    Obj = TpetraMap(num_global, comm); TEST_IERR()
    TEST_ASSERT(Obj%isContiguous())
    call Obj%release(); TEST_IERR()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    Obj = TpetraMap(num_global, elements, comm); TEST_IERR()
    TEST_ASSERT((.not. Obj%isContiguous()))
    call Obj%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isContiguous)

! -------------------------------isDistributed-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isDistributed)
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()

    Obj = TpetraMap(num_global, comm); TEST_IERR()
    if (comm%getSize() == 1) then
      TEST_ASSERT((.not. Obj%isDistributed()))
    else
      TEST_ASSERT(Obj%isDistributed())
    end if
    call Obj%release(); TEST_IERR()

    num_global = 4
    Obj = TpetraMap(num_global, comm, TpetraLocallyReplicated); TEST_IERR()
    TEST_ASSERT((.not. Obj%isDistributed()))
    call Obj%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isDistributed)

! --------------------------------isCompatible-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isCompatible)
    type(TpetraMap) :: Obj1, Obj2
    integer(global_ordinal_type) :: num_global, elements(4)
    integer :: k
    num_global = 4*comm%getSize()

    Obj1 = TpetraMap(num_global, comm); TEST_IERR()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    Obj2 = TpetraMap(num_global, elements, comm); TEST_IERR()

    ! Cyclic map should be compatible
    TEST_ASSERT(Obj1%isCompatible(Obj2))

    call Obj2%release(); TEST_IERR()

    num_global = num_global + 10
    Obj2 = TpetraMap(num_global, comm); TEST_IERR()

    TEST_ASSERT((.not. Obj1%isCompatible(Obj2)))

    call Obj2%release(); TEST_IERR()
    call Obj1%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isCompatible)

! ----------------------------------isSameAs---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isSameAs)
    type(TpetraMap) :: Obj1, Obj2
    integer :: num_local
    integer(global_ordinal_type) :: num_global, elements(4)
    integer :: k
    num_global = 4*comm%getSize()

    Obj1 = TpetraMap(num_global, comm); TEST_IERR()

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      Obj2 = TpetraMap(num_global, elements, comm); TEST_IERR()

      ! Cyclic map should not be SameAs
      TEST_ASSERT((.not. Obj1%isSameAs(Obj2)))

      call Obj2%release(); TEST_IERR()

    end if

    num_local = 4
    Obj2 = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()

    TEST_ASSERT(Obj1%isSameAs(Obj2))

    call Obj1%release(); TEST_IERR()
    call Obj2%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isSameAs)

! -------------------------------locallySameAs-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_locallySameAs)
    type(TpetraMap) :: Obj1, Obj2
    integer :: num_local
    integer(global_ordinal_type) :: num_global, elements(4)
    integer :: k
    num_global = 4*comm%getSize()

    Obj1 = TpetraMap(num_global, comm); TEST_IERR()

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      Obj2 = TpetraMap(num_global, elements, comm); TEST_IERR()

      ! Cyclic map should not be locallySameAs
      TEST_ASSERT((.not. Obj1%locallySameAs(Obj2)))

      call Obj2%release(); TEST_IERR()

    end if

    num_local = 4
    Obj2 = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
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

    Obj = TpetraMap(num_global, comm); TEST_IERR()

    ! We only test the comm returned has the same rank and size.  More
    ! comprehensive testing is (assumed to be) done in Teuchos itself.
    tcomm = Obj%getComm(); TEST_IERR()
    TEST_ASSERT(tcomm%getRank() == comm%getRank())
    TEST_ASSERT(tcomm%getSize() == comm%getSize())

    call Obj%release(); TEST_IERR()
    call tcomm%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getComm)

! --------------------------------description--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_description)
    type(TpetraMap) :: Obj
    character(kind=C_CHAR, len=:), allocatable :: description
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    Obj = TpetraMap(num_global, comm); TEST_IERR()
    description = Obj%description(); TEST_IERR()
    call Obj%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_description)

! ------------------------------isLocallyFitted------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isLocallyFitted)
    type(TpetraMap) :: map1, map2
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    map1 = TpetraMap(num_global, comm); TEST_IERR()
    map2 = TpetraMap(num_global, comm); TEST_IERR()
    TEST_ASSERT(map1%isLocallyFitted(map2))
    call map1%release(); TEST_IERR()
    call map2%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMap_isLocallyFitted)

end program test_TpetraMap
