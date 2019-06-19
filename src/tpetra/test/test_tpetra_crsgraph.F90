! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program test_TpetraCrsGraph
#include "ForTrilinosTpetra_config.hpp"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  implicit none
  type(TeuchosComm) :: comm
  character(len=26), parameter :: FILENAME="test_tpetra_crsgraph.F90"

  SETUP_TEST()

#ifdef HAVE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm()
#endif

  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_ActiveFillLocal)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_Parameters)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_WithColmap)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_TwoArraysESFC)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_SortingTests)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_EmptyFillComplete)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_GetEntities)

!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getComm)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getDomainMap)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getRangeMap)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getImporter)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getExporter)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isLowerTriangular)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isUpperTriangular)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isLocallyIndexed)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isGloballyIndexed)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isStorageOptimized)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_supportsRowViews)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_description)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_replaceColMap)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_replaceDomainMapAndImporter)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_removeEmptyProcessesInPlace)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_haveGlobalConstants)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_computeGlobalConstants)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getGlobalRowCopy)
!  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getGlobalRowView)
!
  call comm%release()

  TEARDOWN_TEST()

contains


  ! ----------------------------- ActiveFill --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_ActiveFillLocal)
    type(TpetraMap) :: Map
    type(ParameterList) :: params
    type(TpetraCrsGraph) :: Graph
    integer :: row, indices(1)

    OUT0("Starting TpetraCrsGraph_ActiveFillLocal!")

    ! create Map
    Map = TpetraMap(TPETRA_GLOBAL_INVALID, 1, comm); TEST_IERR()

    params = ParameterList("ANONYMOUS")
    Graph = TpetraCrsGraph(Map, Map, 1_size_type, TpetraStaticProfile)
    TEST_ASSERT(Graph%isFillActive())
    TEST_ASSERT((.not. Graph%isFillComplete()))

    row = 1
    indices(1) = 1
    call Graph%insertLocalIndices(row, indices)

    call params%set("Optimize Storage", .false.)
    call graph%fillComplete(params)
    TEST_ASSERT((.not. graph%isFillActive()))

    TEST_ASSERT(graph%isFillComplete())
    TEST_THROW(call graph%insertLocalIndices(row, indices))
    TEST_THROW(call graph%removeLocalIndices(row))
    TEST_THROW(call Graph%globalAssemble())
    TEST_THROW(call Graph%fillComplete())

    call Graph%resumeFill()
    TEST_ASSERT(graph%isFillActive())
    TEST_ASSERT((.not. graph%isFillComplete()))
    TEST_ASSERT(graph%getProfileType() == TpetraStaticProfile)
    TEST_NOTHROW(call graph%insertLocalIndices(row, indices))

    TEST_NOTHROW(call graph%fillComplete())
    TEST_ASSERT((.not. graph%isFillActive()))
    TEST_ASSERT(graph%isFillComplete())

    call Graph%release()
    call Map%release()
    call params%release()

    OUT0("Finished TpetraCrsGraph_ActiveFillLocal!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_ActiveFillLocal)

  ! ----------------------------getValidParameters---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_Parameters)
    type(TpetraMap) :: Map
    type(ParameterList) :: params
    type(TpetraCrsGraph) :: Graph
    integer :: row, indices(1)

    OUT0("Starting TpetraCrsGraph_Parameters!")

    ! create Map
    Map = TpetraMap(TPETRA_GLOBAL_INVALID, 1, comm); TEST_IERR()

    params = ParameterList("ANONYMOUS")
    Graph = TpetraCrsGraph(Map, Map, 1_size_type, TpetraStaticProfile)
    params = Graph%getValidParameters()
    call Graph%setParameterList(params)

    row = 1
    indices(1) = 1
    call Graph%insertLocalIndices(row, indices)
    call graph%fillComplete(params)

    call params%release()
    call Graph%release()
    call Map%release()

    OUT0("Finished TpetraCrsGraph_Parameters!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_Parameters)

  ! ------------------------------- WithColMap ------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_WithColMap)
    type(TpetraMap) :: rmap, cmap, tmpmap
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: num_procs, num_ent_per_row
    integer :: num_local, lclrow
    integer(global_ordinal_type) :: myrowind
    integer(global_ordinal_type), allocatable :: indices(:)

    OUT0("Starting TpetraCrsGraph_WithColMap!")

    num_procs = comm%getSize()
    if (num_procs == 1) return

    num_local = 1
    rmap = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    cmap = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    ! must allocate enough for all submitted indices.
    num_ent_per_row = 2
    Graph = TpetraCrsGraph(rmap, cmap, num_ent_per_row, TpetraStaticProfile)
    TEST_ASSERT(Graph%hasColMap())
    lclrow = 1
    myrowind = rmap%getGlobalElement(lclrow);

    ! insertGlobalIndices doesn't do column Map filtering, so we have to test
    ! whether each of the column indices to insert is invalid.
    if (cmap%isNodeGlobalElement(myrowind)) then
      if (cmap%isNodeGlobalElement(myrowind+1)) then
        allocate(indices(2))
        indices = [myrowind, myrowind+1]
        TEST_NOTHROW(call Graph%insertGlobalIndices(myrowind, indices))
        deallocate(indices)
      else
        allocate(indices(1))
        indices(1) = myrowind
        TEST_NOTHROW(call Graph%insertGlobalIndices(myrowind, indices))
        deallocate(indices)
      end if
    end if
    TEST_NOTHROW(call Graph%fillComplete())

    tmpmap = TpetraMap()
    tmpmap = Graph%getRowMap()
    TEST_EQUALITY(tmpmap%getNodeNumElements(), rmap%getNodeNumElements())
    tmpmap = Graph%getColMap()
    TEST_EQUALITY(tmpmap%getNodeNumElements(), cmap%getNodeNumElements())
    TEST_EQUALITY(Graph%getNumEntriesInGlobalRow(myrowind), 1_size_type)

    call tmpmap%release()
    call Graph%release()
    call rmap%release()
    call cmap%release()

    OUT0("Finished TpetraCrsGraph_WithColMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_WithColMap)

  ! ------------------------------- TwoArraysESFC ---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_TwoArraysESFC)
    type(TpetraMap) :: rmap, tmpmap
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: num_procs
    integer :: num_local
    integer(size_type), allocatable :: rowptr(:)
    integer, allocatable :: colind(:)

    OUT0("Starting TpetraCrsGraph_TwoArraysESFC!")

    num_procs = comm%getSize()
    if (num_procs == 1) return

    num_local = 2;
    rmap = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()

    allocate(rowptr(num_local+1))
    allocate(colind(num_local)) ! one unknown per row
    rowptr(1:3) = [1, 2, 3]
    colind(1:2) = [1, 2]

    Graph = TpetraCrsGraph(rmap, rmap, rowptr, colind)
    TEST_ASSERT(Graph%hasColMap())

    TEST_NOTHROW(call Graph%expertStaticFillComplete(rmap, rmap))

    tmpmap = TpetraMap()
    tmpmap = Graph%getRowMap()
    TEST_EQUALITY(tmpmap%getNodeNumElements(), rmap%getNodeNumElements())

    tmpmap = Graph%getColMap()
    TEST_EQUALITY(tmpmap%getNodeNumElements(), rmap%getNodeNumElements())

    deallocate(rowptr, colind)
    call tmpmap%release()
    call Graph%release()
    call rmap%release()

    OUT0("Finished TpetraCrsGraph_TwoArraysESFC!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_TwoArraysESFC)

#if 0
  ! TODO: Implement setAllIndices
  ! ------------------------------ SetAllIndices ----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_SetAllIndices)
    type(TpetraMap) :: rmap, tmpmap
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: num_procs
    integer :: num_local
    integer(size_type), allocatable :: rowptr(:)
    integer, allocatable :: colind(:)
    integer(global_ordinal_type), allocatable :: indices(:)

    OUT0("Starting TpetraCrsGraph_SetAllIndices!")
    num_procs = comm%getSize()
    !if (num_procs <= 1) return

    num_local = 2
    rmap = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm)
    allocate(rowptr(num_local+1))
    allocate(colind(num_local)) ! one unknown per row
    rowptr(1:3) = [1, 2, 3]
    colind(1:2) = [1, 2]

    Graph = TpetraCrsGraph(rmap, rmap, 0_size_type, TpetraStaticProfile)

    TEST_NOTHROW(call Graph%setAllIndices(rowptr, colind))
    TEST_ASSERT(Graph%hasColMap())

    TEST_NOTHROW(call Graph%expertStaticFillComplete(rmap,rmap))

    tmpmap = TpetraMap()
    tmpmap = Graph%getRowMap()
    TEST_EQUALITY(tmpmap%getNodeNumElements(), rmap%getNodeNumElements())

    tmpmap = Graph%getColMap()
    TEST_EQUALITY(tmpmap%getNodeNumElements(), rmap%getNodeNumElements())

    deallocate(rowptr, colind)
    call tmpmap%release()
    call Graph%release()
    call rmap%release()

    OUT0("Finished TpetraCrsGraph_SetAllIndices!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_SetAllIndices)
#endif

  ! ------------------------------ SortingTests ------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_SortingTests)
    type(TpetraMap) :: map
    type(TpetraCrsGraph) :: Graph
    type(ParameterList) :: params
    integer(size_type) :: kk, nument
    integer :: num_local
    logical :: sorting_check
    integer(global_ordinal_type) :: j, jj, jstart, jfinish, firstind, lastind
    integer(global_ordinal_type), allocatable :: jinds(:)
    integer :: k, kstart, kfinish
    integer, allocatable :: kinds(:)

    OUT0("Starting TpetraCrsGraph_SortingTests!")

    ! create a Map
    num_local = 10;
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm)

    nument = 4
    Graph = TpetraCrsGraph(map, map, nument*num_local); TEST_IERR()
    TEST_ASSERT(graph%isSorted())

    ! insert entries; shouldn't be sorted anymore
    jstart = map%getMinGlobalIndex()
    jfinish = map%getMaxGlobalIndex()
    jj = 1
    allocate(jinds(100))
    do j = jstart, jfinish
      firstind = mod(j+5, map%getMaxAllGlobalIndex())
      if (map%isNodeGlobalElement(firstind)) then
        jinds(jj) = firstind; jj = jj + 1
      end if
      jinds(jj) = j; jj = jj + 1
      lastind = mod((j-5), map%getMaxAllGlobalIndex())
      if (map%isNodeGlobalElement(lastind)) then
        jinds(jj) = lastind; jj = jj + 1
      end if
      call Graph%insertGlobalIndices(j, jinds(1:jj-1)); TEST_IERR()
    end do
    TEST_ASSERT((.not. Graph%isSorted()))
    deallocate(jinds)

    ! fill complete; should be sorted now
    params = ParameterList("ANONYMOUS")
    call params%set("Optimize Storage", .false.)
    call graph%fillComplete(params)

    sorting_check = .true.
    kstart = map%getMinLocalIndex()
    kfinish = map%getMaxLocalIndex()
    do k = kstart, kfinish
      nument = Graph%getNumEntriesInLocalRow(k)
      allocate(kinds(nument))
      write (*,*) "Allocated", nument, "for k=", k
      call Graph%getLocalRowCopy(k, kinds, nument); TEST_IERR()
      do kk = 2, nument
        if (kinds(kk-1) > kinds(kk)) then
          sorting_check = .false.
          exit
        endif
      end do
      TEST_ASSERT((sorting_check .eqv. graph%isSorted()))
      deallocate(kinds)
    end do

    ! resume fill; should still be sorted
    call Graph%resumeFill(); TEST_IERR()
    TEST_ASSERT(graph%isSorted())

    sorting_check = .true.
    kstart = map%getMinLocalIndex()
    kfinish = map%getMaxLocalIndex()
    do k = kstart, kfinish
      nument = Graph%getNumEntriesInLocalRow(k)
      allocate(kinds(nument))
      call Graph%getLocalRowCopy(k, kinds, nument); TEST_IERR()
      do kk = 2, nument
        if (kinds(kk-1) > kinds(kk)) then
          sorting_check = .false.
          exit
        endif
      end do
      TEST_ASSERT((sorting_check .eqv. graph%isSorted()))
      deallocate(kinds)
    end do

    ! insert a column-index; currently, this invalidates sorting, though it may
    ! change in the future
    k = 1
    allocate(kinds(1))
    kinds(1) = 1
    call Graph%insertLocalIndices(k, kinds); TEST_IERR()
    TEST_ASSERT((.not. graph%isSorted()))
    deallocate(kinds)

    ! fill complete, check one more time
    call params%set("Optimize Storage", .false.)
    call graph%fillComplete(params)

    sorting_check = .true.
    kstart = map%getMinLocalIndex()
    kfinish = map%getMaxLocalIndex()
    do k = kstart, kfinish
      nument = Graph%getNumEntriesInLocalRow(k)
      allocate(kinds(nument))
      call Graph%getLocalRowCopy(k, kinds, nument); TEST_IERR()
      do kk = 2, nument
        if (kinds(kk-1) > kinds(kk)) then
          sorting_check = .false.
          exit
        endif
      end do
      TEST_ASSERT((sorting_check .eqv. graph%isSorted()))
      deallocate(kinds)
    end do

    call Graph%release()
    call map%release()
    call params%release()
    OUT0("Finished TpetraCrsGraph_SortingTests!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_SortingTests)

  ! --------------------------- EmptyFillComplete ---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_EmptyFillComplete)
    type(TpetraMap) :: map
    type(TpetraCrsGraph) :: Graph
    integer :: num_local

    OUT0("Starting TpetraCrsGraph_EmptyFillComplete!")

    ! create a Map
    num_local = 10
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm)

    ! create static-profile graph, fill-complete without inserting
    ! (and therefore, without allocating)
    Graph = TpetraCrsGraph(map, 1_size_type, TpetraStaticProfile)
    call Graph%fillComplete()
    call Graph%release()

    call map%release()

    OUT0("Finished TpetraCrsGraph_EmptyFillComplete!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_EmptyFillComplete)

  ! ------------------------------- GetEntities ------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_GetEntities)
    type(TpetraMap) :: rmap, cmap
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: num_procs, num_ent_per_row, i_LO
    integer :: num_local, lclrow
    integer(global_ordinal_type) :: myrowind
    integer(global_ordinal_type), allocatable :: indices(:)
    integer :: commsize

    OUT0("Starting TpetraCrsGraph_GetEntities!")

    num_procs = comm%getSize()
    if (num_procs == 1) return

    num_local = 1
    rmap = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    cmap = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    ! must allocate enough for all submitted indices.
    num_ent_per_row = 2
    Graph = TpetraCrsGraph(rmap, cmap, num_ent_per_row, TpetraStaticProfile)
    TEST_ASSERT(Graph%hasColMap())
    lclrow = 1
    myrowind = rmap%getGlobalElement(lclrow);

    ! insertGlobalIndices doesn't do column Map filtering, so we have to test
    ! whether each of the column indices to insert is invalid.
    if (cmap%isNodeGlobalElement(myrowind)) then
      if (cmap%isNodeGlobalElement(myrowind+1)) then
        allocate(indices(2))
        indices = [myrowind, myrowind+1]
        TEST_NOTHROW(call Graph%insertGlobalIndices(myrowind, indices))
        deallocate(indices)
      else
        allocate(indices(1))
        indices(1) = myrowind
        TEST_NOTHROW(call Graph%insertGlobalIndices(myrowind, indices))
        deallocate(indices)
      end if
    end if
    TEST_NOTHROW(call Graph%fillComplete())

    commsize = comm%getSize()
    i_LO = Graph%getNodeNumRows()
    TEST_EQUALITY(i_LO, 1_size_type)
    i_LO = Graph%getNodeNumCols()
    TEST_EQUALITY(i_LO, 1_size_type)
    i_LO = Graph%getGlobalNumEntries()
    TEST_EQUALITY(i_LO, int(commsize, kind=size_type))
    i_LO = Graph%getNodeNumEntries()
    TEST_EQUALITY(i_LO, 1_size_type)
    i_LO = Graph%getNodeAllocationSize()
    TEST_EQUALITY(i_LO, 1_size_type)
    i_LO = Graph%getNumAllocatedEntriesInGlobalRow(myrowind)
    TEST_EQUALITY(i_LO, 1_size_type)
    i_LO = Graph%getNumAllocatedEntriesInLocalRow(lclrow)
    TEST_EQUALITY(i_LO, 1_size_type)
    i_LO = Graph%getGlobalMaxNumRowEntries()
    TEST_EQUALITY(i_LO, 1_size_type)
    i_LO = Graph%getNodeMaxNumRowEntries()
    TEST_EQUALITY(i_LO, 1_size_type)

    call Graph%release()
    call rmap%release()
    call cmap%release()

    OUT0("Finished TpetraCrsGraph_GetEntities!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_GetEntities)

#if 0
  ! ---------------------------------getComm---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getComm)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_getComm!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%getComm(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_getComm: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_getComm!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getComm)

  ! -------------------------------getDomainMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getDomainMap)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_getDomainMap!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%getDomainMap(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_getDomainMap: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_getDomainMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getDomainMap)

  ! -------------------------------getRangeMap-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getRangeMap)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_getRangeMap!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%getRangeMap(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_getRangeMap: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_getRangeMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getRangeMap)

  ! -------------------------------getImporter-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getImporter)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_getImporter!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%getImporter(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_getImporter: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_getImporter!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getImporter)

  ! -------------------------------getExporter-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getExporter)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_getExporter!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%getExporter(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_getExporter: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_getExporter!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getExporter)

  ! ----------------------------isLowerTriangular----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isLowerTriangular)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_isLowerTriangular!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%isLowerTriangular(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_isLowerTriangular: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_isLowerTriangular!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isLowerTriangular)

  ! ----------------------------isUpperTriangular----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isUpperTriangular)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_isUpperTriangular!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%isUpperTriangular(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_isUpperTriangular: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_isUpperTriangular!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isUpperTriangular)

  ! -----------------------------isLocallyIndexed----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isLocallyIndexed)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_isLocallyIndexed!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%isLocallyIndexed(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_isLocallyIndexed: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_isLocallyIndexed!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isLocallyIndexed)

  ! ----------------------------isGloballyIndexed----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isGloballyIndexed)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_isGloballyIndexed!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%isGloballyIndexed(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_isGloballyIndexed: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_isGloballyIndexed!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isGloballyIndexed)

  ! ----------------------------isStorageOptimized---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isStorageOptimized)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_isStorageOptimized!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%isStorageOptimized(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_isStorageOptimized: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_isStorageOptimized!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isStorageOptimized)

  ! -----------------------------supportsRowViews----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_supportsRowViews)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_supportsRowViews!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%supportsRowViews(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_supportsRowViews: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_supportsRowViews!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_supportsRowViews)

  ! -------------------------------description-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_description)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_description!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%description(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_description: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_description!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_description)

  ! ------------------------------replaceColMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_replaceColMap)
    type(TpetraCrsGraph) :: Obj
    type(TpetraMap) :: newcolmap
    OUT0("Starting TpetraCrsGraph_replaceColMap!")

    success = .false.

    !newcolmap = TpetraMap(); TEST_IERR()
    !Obj = TpetraCrsGraph(); TEST_IERR()
    !call Obj%replaceColMap(newcolmap); TEST_IERR()

    !call newcolmap%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_replaceColMap: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_replaceColMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_replaceColMap)

  ! -----------------------replaceDomainMapAndImporter------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_replaceDomainMapAndImporter)
    type(TpetraCrsGraph) :: Obj
    type(TpetraMap) :: newdomainmap
    type(TpetraImport) :: newimporter
    OUT0("Starting TpetraCrsGraph_replaceDomainMapAndImporter!")

    success = .false.

    !newdomainmap = TpetraMap(); TEST_IERR()
    !newimporter = TpetraImport(); TEST_IERR()
    !Obj = TpetraCrsGraph(); TEST_IERR()
    !call Obj%replaceDomainMapAndImporter(newdomainmap, newimporter); TEST_IERR()

    !call newdomainmap%release(); TEST_IERR()
    !call newimporter%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_replaceDomainMapAndImporter: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_replaceDomainMapAndImporter!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_replaceDomainMapAndImporter)

  ! -----------------------removeEmptyProcessesInPlace------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_removeEmptyProcessesInPlace)
    type(TpetraCrsGraph) :: Obj
    type(TpetraMap) :: newmap
    OUT0("Starting TpetraCrsGraph_removeEmptyProcessesInPlace!")

    success = .false.

    !newmap = TpetraMap(); TEST_IERR()
    !Obj = TpetraCrsGraph(); TEST_IERR()
    !call Obj%removeEmptyProcessesInPlace(newmap); TEST_IERR()

    !call newmap%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_removeEmptyProcessesInPlace: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_removeEmptyProcessesInPlace!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_removeEmptyProcessesInPlace)

  ! ---------------------------haveGlobalConstants---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_haveGlobalConstants)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_haveGlobalConstants!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !fresult = Obj%haveGlobalConstants(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_haveGlobalConstants: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_haveGlobalConstants!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_haveGlobalConstants)

  ! --------------------------computeGlobalConstants-------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_computeGlobalConstants)
    type(TpetraCrsGraph) :: Obj
    OUT0("Starting TpetraCrsGraph_computeGlobalConstants!")

    success = .false.

    !Obj = TpetraCrsGraph(); TEST_IERR()
    !call Obj%computeGlobalConstants(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_computeGlobalConstants: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_computeGlobalConstants!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_computeGlobalConstants)

  ! ----------------------------insertLocalIndices---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_insertLocalIndices)
    type(TpetraCrsGraph) :: Obj
    integer(C_INT) :: localrow
    integer(C_INT), allocatable :: indices(:)
    OUT0("Starting TpetraCrsGraph_insertLocalIndices!")

    success = .false.

    localrow = 0
    !allocate(indices(:)(0))
    !Obj = TpetraCrsGraph(); TEST_IERR()
    !call Obj%insertLocalIndices(localrow, indices(:)); TEST_IERR()

    !deallocate(indices(:))
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_insertLocalIndices: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_insertLocalIndices!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_insertLocalIndices)

  ! -----------------------------getGlobalRowCopy----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalRowCopy)
    type(TpetraCrsGraph) :: Obj
    integer(C_LONG_LONG) :: globalrow
    integer(C_LONG_LONG), allocatable :: indices(:)
    integer(C_SIZE_T) :: numindices
    OUT0("Starting TpetraCrsGraph_getGlobalRowCopy!")

    success = .false.

    globalrow = 0
    !allocate(indices(:)(0))
    numindices = 0
    !Obj = TpetraCrsGraph(); TEST_IERR()
    !call Obj%getGlobalRowCopy(globalrow, indices(:), numindices); TEST_IERR()

    !deallocate(indices(:))
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_getGlobalRowCopy: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_getGlobalRowCopy!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalRowCopy)

  ! ------------------------------getGlobalRowView------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalRowView)
    type(TpetraCrsGraph) :: Obj
    integer(C_LONG_LONG) :: gblrow
    integer(C_LONG_LONG), allocatable :: lclcolinds(:)
    OUT0("Starting TpetraCrsGraph_getGlobalRowView!")

    success = .false.

    gblrow = 0
    !allocate(lclcolinds(:)(0))
    !Obj = TpetraCrsGraph(); TEST_IERR()
    !call Obj%getGlobalRowView(gblrow, lclcolinds(:)); TEST_IERR()

    !deallocate(lclcolinds(:))
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsGraph_getGlobalRowView: Test not yet implemented'

    OUT0("Finished TpetraCrsGraph_getGlobalRowView!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalRowView)
#endif

end program test_TpetraCrsGraph
