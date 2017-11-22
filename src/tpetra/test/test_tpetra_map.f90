program test_TpetraMap
#include "ForTrilinosTpetra_config.hpp"
#include "FortranTestMacros.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  DECLARE_TEST_VARIABLES()
  type(TeuchosComm) :: comm
  integer(global_size_type), parameter :: invalid=-1
  integer(global_ordinal_type), parameter :: index_base=1

  INITIALIZE_TEST()

#ifdef HAVE_MPI
  call comm%create(MPI_COMM_WORLD)
  CHECK_IERR()
#else
  call comm%create()
  CHECK_IERR()
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
  CHECK_IERR()

  SHUTDOWN_TEST()

contains

! ---------------------------------isOneToOne--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isOneToOne)
    integer :: jerr
    logical(c_bool) :: bool
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global, indices(4)
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    bool = Obj%isOneToOne(); TEST_IERR()
    if (.not. bool) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isOneToOne: Expected map to be one to one"
    end if

    call Obj%release(); TEST_IERR()

    if (comm%getSize() > 1) then
      indices = [1, 2, 3, 4]
      call Obj%create(num_global, indices, index_base, comm); TEST_IERR()

      bool = Obj%isOneToOne(); TEST_IERR()
      if (bool) then
        jerr = jerr + 1
        if (comm%getRank() == 0) &
          write(*,*) "isOneToOne: Expected map to NOT be one to one"
      end if

      call Obj%release(); TEST_IERR()

    end if

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isOneToOne)

! ----------------------------getGlobalNumElements---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getGlobalNumElements)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%getGlobalNumElements(); TEST_IERR()
    if (fresult /= num_global) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getGlobalNumElements: Expected ", num_global, " elements, got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getGlobalNumElements)

! -----------------------------getNodeNumElements----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getNodeNumElements)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_SIZE_T) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%getNodeNumElements()
    if (fresult /= 4) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getNodeNumElements: Expected ", 4, " elements, got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getNodeNumElements)

! --------------------------------getIndexBase-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getIndexBase)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base_0
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%getIndexBase(); TEST_IERR()
    if (fresult /= 1) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getIndexBase: Expected indexBase = ", 1, " got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    index_base_0 = 0
    call Obj%create(num_global, index_base_0, comm); TEST_IERR()

    fresult = Obj%getIndexBase(); TEST_IERR()
    if (fresult /= 0) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getIndexBase: Expected indexBase = ", 0, " got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getIndexBase)

! ------------------------------getMinLocalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMinLocalIndex)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%getMinLocalIndex(); TEST_IERR()
    if (fresult /= 1) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMinLocalIndex: Expected minLocalIndex = ", 1, " got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMinLocalIndex)

! ------------------------------getMaxLocalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMaxLocalIndex)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%getMaxLocalIndex(); TEST_IERR()
    if (fresult /= 4) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMaxLocalIndex: Expected maxLocalIndex = ", 4, " got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMaxLocalIndex)

! -----------------------------getMinGlobalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMinGlobalIndex)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%getMinGlobalIndex(); TEST_IERR()
    expected = comm%getRank() * 4 + 1
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMindGlobalIndex: Expected minGlobalIndex = ", expected, " got ", fresult
    end if
    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMinGlobalIndex)

! -----------------------------getMaxGlobalIndex------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMaxGlobalIndex)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%getMaxGlobalIndex(); TEST_IERR()
    expected = comm%getRank() * 4 + 4
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMaxGlobalIndex: Expected maxGloblIndex = ", expected, " got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMaxGlobalIndex)

! ----------------------------getMinAllGlobalIndex---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMinAllGlobalIndex)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%getMinAllGlobalIndex(); TEST_IERR()
    if (fresult /= index_base) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMinAllGlobalIndex: Expected minAllGlobalIndex = ", index_base, " got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMinAllGlobalIndex)

! ----------------------------getMaxAllGlobalIndex---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getMaxAllGlobalIndex)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%getMaxAllGlobalIndex(); TEST_IERR()
    if (fresult /= num_global) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMaxAllGlobalIndex: Expected maxAllGlobalIndex = ", num_global, " got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getMaxAllGlobalIndex)

! ------------------------------getLocalElement------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getLocalElement)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    globalindex = comm%getRank() * 4 + 1
    fresult = Obj%getLocalElement(globalindex); TEST_IERR()
    if (fresult /= 1) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getLocalElement: Expected local element = ", 1, " got ", fresult
    end if

    globalindex = comm%getRank() * 4 + 4
    fresult = Obj%getLocalElement(globalindex); TEST_IERR()
    if (fresult /= 4) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getLocalElement: Expected local element = ", 4, " got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getLocalElement)

! ------------------------------getGlobalElement------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMap_getGlobalElement)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    localindex = 1
    expected = comm%getRank() * 4 + 1
    fresult = Obj%getGlobalElement(localindex); TEST_IERR()
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getGlobalElement: Expected global element = ", expected, " got ", fresult
    end if

    localindex = 4
    expected = comm%getRank() * 4 + 4
    fresult = Obj%getGlobalElement(localindex); TEST_IERR()
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getGlobalElement: Expected global element = ", expected, " got ", fresult
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

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
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    localindex = 1
    fresult = Obj%isNodeLocalElement(localindex); TEST_IERR()
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeLocalElement: expected 1 to be a local index"
    end if

    localindex = 5
    fresult = Obj%isNodeLocalElement(localindex); TEST_IERR()
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeLocalElement: expected 5 to NOT be a local index"
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isNodeLocalElement)

! ----------------------------isNodeGlobalElement----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isNodeGlobalElement)
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    globalindex = int(comm%getRank() * 4 + 1, kind=global_ordinal_type)
    fresult = Obj%isNodeGlobalElement(globalindex); TEST_IERR()
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeGlobalElement: expected 1 to be a global index"
    end if
    globalindex = num_global + 1

    fresult = Obj%isNodeGlobalElement(globalindex); TEST_IERR()
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeGlobalElement: expected 5 to NOT be a global index"
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isNodeGlobalElement)

! ---------------------------------isUniform---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isUniform)
    integer :: jerr
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer :: k
    integer(global_ordinal_type) :: num_global, elements(4)
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%isUniform(); TEST_IERR()
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isUniform: expected map to be uniform"
    end if

    call Obj%release(); TEST_IERR()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm); TEST_IERR()

    fresult = Obj%isUniform(); TEST_IERR()
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isUniform: expected map to NOT be uniform"
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isUniform)

! --------------------------------isContiguous-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isContiguous)
    integer :: jerr
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer :: k
    integer(global_ordinal_type) :: num_global, elements(4)
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%isContiguous(); TEST_IERR()
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isContiguous: expected map to be uniform"
    end if

    call Obj%release(); TEST_IERR()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm); TEST_IERR()

    fresult = Obj%isContiguous(); TEST_IERR()
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isContiguous: expected map to NOT be uniform"
    end if

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isContiguous)

! -------------------------------isDistributed-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isDistributed)
    integer :: jerr
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, elements(4)
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj%isDistributed(); TEST_IERR()
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

    call Obj%release(); TEST_IERR()

    ! All elements have 4 entries and map has only 4 entries so the map should
    ! not be distributed
    num_global = 4
    elements = [1, 2, 3, 4]
    call Obj%create(num_global, elements, index_base, comm); TEST_IERR()

    fresult = Obj%isDistributed(); TEST_IERR()
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

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isDistributed)

! --------------------------------isCompatible-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isCompatible)
    integer :: jerr
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, elements(4)
    integer :: k
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj1%create(num_global, index_base, comm); TEST_IERR()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj2%create(num_global, elements, index_base, comm); TEST_IERR()

    ! Cyclic map should be compatible
    fresult = Obj1%isCompatible(Obj2); TEST_IERR()
    if (.not. fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isCompatible: Expected maps to be compatible"
    end if

    call Obj2%release(); TEST_IERR()

    num_global = num_global + 10
    call Obj2%create(num_global, index_base, comm); TEST_IERR()

    fresult = Obj1%isCompatible(Obj2); TEST_IERR()
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "isCompatible: Expected maps to NOT be compatible"
    end if

    call Obj2%release(); TEST_IERR()
    call Obj1%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isCompatible)

! ----------------------------------isSameAs---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_isSameAs)
    integer :: jerr
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    logical(C_BOOL) :: fresult
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, elements(4)
    integer :: k
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj1%create(num_global, index_base, comm); TEST_IERR()

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      call Obj2%create(num_global, elements, index_base, comm); TEST_IERR()

      ! Cyclic map should not be SameAs
      fresult = Obj1%isSameAs(Obj2); TEST_IERR()
      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isSameAs: Expected maps to NOT be same"
      end if

      call Obj2%release(); TEST_IERR()

    end if

    num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm); TEST_IERR()

    fresult = Obj1%isSameAs(Obj2); TEST_IERR()
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "isSameAs: Expected maps to be same"
    end if

    call Obj1%release(); TEST_IERR()
    call Obj2%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_isSameAs)

! -------------------------------locallySameAs-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_locallySameAs)
    integer :: jerr
    type(TpetraMap) :: Obj1, Obj2
    logical(C_BOOL) :: fresult
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, elements(4)
    integer :: k
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj1%create(num_global, index_base, comm); TEST_IERR()

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      call Obj2%create(num_global, elements, index_base, comm); TEST_IERR()

      ! Cyclic map should not be locallySameAs
      fresult = Obj1%locallySameAs(Obj2); TEST_IERR()
      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "locallySameAs: Expected maps to NOT be same"
      end if

      call Obj2%release(); TEST_IERR()

    end if

    num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm); TEST_IERR()

    fresult = Obj1%locallySameAs(Obj2); TEST_IERR()
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "locallySameAs: Expected maps to be same"
    end if

    call Obj1%release(); TEST_IERR()
    call Obj2%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_locallySameAs)

! ----------------------------------getComm----------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_getComm)
    integer :: jerr
    type(TpetraMap) :: Obj
    type(TeuchosComm) :: fresult
    integer(global_ordinal_type) :: num_global
    jerr = 0
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm); TEST_IERR()

    ! We only test the comm returned has the same rank and size.  More
    ! comprehensive testing is (assumed to be) done in Teuchos itself.
    fresult = Obj%getComm(); TEST_IERR()
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

    call Obj%release(); TEST_IERR()

    success = (jerr == 0)

  END_FORTRILINOS_UNIT_TEST(TpetraMap_getComm)

! --------------------------------description--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMap_description)
    type(TpetraMap) :: Obj
    type(string) :: fresult
    integer(global_ordinal_type) :: num_global
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm); TEST_IERR()
    fresult = Obj%description(); TEST_IERR()
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
