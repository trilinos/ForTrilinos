! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program test_TpetraCrsGraph
#include "ForTrilinos_config.h"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra
  use test_tpetra_crsgraph_helper

  implicit none
  type(TeuchosComm) :: comm
  character(len=26), parameter :: FILENAME="test_tpetra_crsgraph.F90"

  SETUP_TEST()

#if FORTRILINOS_USE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm()
#endif

  !Fatter tests
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_ActiveFillLocal)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_WithColMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_TwoArraysESFC)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_SortingTests)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_EmptyFillComplete)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_GetEntities)

  !Unit tests consistent with scripts/autogen_tests.py
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isIdenticalTo)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_insertGlobalIndices)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_insertLocalIndices)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_removeLocalIndices)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getComm)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getRowMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getColMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getDomainMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getRangeMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getImporter)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getExporter)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getGlobalNumRows)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getGlobalNumCols)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getGlobalNumEntries)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getGlobalMaxNumRowEntries)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNodeNumRows)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNodeNumCols)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNodeNumEntries)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNodeMaxNumRowEntries)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNumEntriesInGlobalRow)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNumEntriesInLocalRow)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNodeAllocationSize)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNumAllocatedEntriesInGlobalRow)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNumAllocatedEntriesInLocalRow)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getProfileType)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getGlobalRowCopy)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getLocalRowCopy)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNodeRowPtrs)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getNodePackedIndices)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_hasColMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isGloballyIndexed)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isLocallyIndexed)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isFillComplete)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isFillActive)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isSorted)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_isStorageOptimized)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_supportsRowViews)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_description)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_replaceColMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_replaceDomainMapAndImporter)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_removeEmptyProcessesInPlace)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_haveGlobalConstants)
  ADD_SUBTEST_AND_RUN(TpetraCrsGraph_computeGlobalConstants)

  ! In the wrapper but not wanting to really expose for users
  !ADD_SUBTEST_AND_RUN(TpetraCrsGraph_setParameterList)
  !ADD_SUBTEST_AND_RUN(TpetraCrsGraph_getValidParameters)
  !ADD_SUBTEST_AND_RUN(TpetraCrsGraph_globalAssemble)

  ! Deprecated
  call comm%release(); TEST_IERR()

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
    TEST_THROW(call Graph%fillComplete())

    call Graph%resumeFill()
    TEST_ASSERT(graph%isFillActive())
    TEST_ASSERT((.not. graph%isFillComplete()))
    TEST_ASSERT(graph%getProfileType() == TpetraStaticProfile)
    TEST_NOTHROW(call graph%insertLocalIndices(row, indices))

    TEST_NOTHROW(call graph%fillComplete())
    TEST_ASSERT((.not. graph%isFillActive()))
    TEST_ASSERT(graph%isFillComplete())

    call Graph%release(); TEST_IERR()
    call Map%release(); TEST_IERR()
    call params%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_ActiveFillLocal!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_ActiveFillLocal)

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

    call tmpmap%release(); TEST_IERR()
    call Graph%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()
    call cmap%release(); TEST_IERR()

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
    call tmpmap%release(); TEST_IERR()
    call Graph%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_TwoArraysESFC!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_TwoArraysESFC)


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
      !write (*,*) "Allocated", nument, "for k=", k
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

    call Graph%release(); TEST_IERR()
    call map%release(); TEST_IERR()
    call params%release(); TEST_IERR()
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
    call Graph%release(); TEST_IERR()

    call map%release(); TEST_IERR()

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

    call Graph%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()
    call cmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_GetEntities!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_GetEntities)

  ! ------------------------------isIdenticalTo------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isIdenticalTo)
    type(TpetraMap) :: Map
    type(TpetraCrsGraph) :: Graph1, Graph2
    integer(global_ordinal_type) :: num_global
    integer(size_type), parameter :: izero=0, ione=1
    logical :: fresult

    OUT0("Starting TpetraCrsGraph_isIdenticalTo!")

    ! create Map and Graph
    num_global = 4*comm%getSize()
    Map = TpetraMap(num_global, comm); TEST_IERR()
    Graph1 = TpetraCrsGraph(Map, Map, ione, TpetraStaticProfile)
    Graph2 = TpetraCrsGraph(Map, Map, ione, TpetraStaticProfile)

    call Graph1%fillComplete()
    call Graph2%fillComplete()

    TEST_ASSERT(Graph1%isIdenticalTo(Graph2))

    call Graph1%release(); TEST_IERR()
    call Graph2%release(); TEST_IERR()
    call Map%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_isIdenticalTo!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isIdenticalTo)

  ! ----------------------------insertLocalIndices---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_insertLocalIndices)
    type(TpetraCrsGraph) :: Graph
    type(TpetraMap) :: rowMap, colMap
    integer(int_type) :: lclNumRows, lclNumCols, lclRow
    integer(global_ordinal_type) gblNumRows
    integer(int_type) :: numProcs
    integer(int_type), dimension(:), allocatable :: lclColInds
    integer(int_type), dimension(:), allocatable :: curLclColInds
    integer(size_type), parameter :: max_entries_per_row=3
    integer(size_type) :: numEnt=0

    OUT0("Starting TpetraCrsGraph_insertLocalIndices!")

    numProcs = comm%getSize()
    lclNumRows = 10
    gblNumRows = lclNumRows * numProcs
    rowMap = TpetraMap(gblNumRows, lclNumRows, comm)
    lclNumCols = lclNumRows
    colMap = rowMap
    Graph = TpetraCrsGraph(rowMap, colMap, max_entries_per_row, TpetraStaticProfile)

    allocate(lclColInds(2))

    do lclRow = 1, lclNumRows
      lclColInds(1) = 1 + modulo(lclRow + 0, lclNumCols)
      lclColInds(2) = 1 + modulo(lclRow + 1, lclNumCols)
      TEST_NOTHROW(call Graph%insertLocalIndices(lclRow, lclColInds))
    end do

    ! Make the arrays bigger than necessary, just to make sure that
    ! the methods behave correctly.
    allocate(curLclColInds(5))

    do lclRow = 1, lclNumRows
      TEST_NOTHROW(call Graph%getLocalRowCopy(lclRow, curLclColInds, numEnt))
      TEST_EQUALITY(int(numEnt), 2)

      TEST_EQUALITY(curLclColInds(1), 1 + modulo(lclRow + 0, lclNumCols));
      TEST_EQUALITY(curLclColInds(2), 1 + modulo(lclRow + 1, lclNumCols));
    enddo

    call Graph%release(); TEST_IERR()
    call rowMap%release(); TEST_IERR()
    call colMap%release(); TEST_IERR()
    deallocate(lclColInds)
    deallocate(curLclColInds)

    OUT0("Finished TpetraCrsGraph_insertLocalIndices!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_insertLocalIndices)

  ! ----------------------------insertGlobalIndices---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_insertGlobalIndices)
    type(TpetraCrsGraph) :: Graph
    type(TpetraMap) :: rmap
    integer :: irow, nnz
    integer(size_type) :: nument;
    integer(global_ordinal_type), allocatable :: cols_expected(:)
    integer(global_ordinal_type), allocatable :: cols(:)
    integer(global_ordinal_type) :: gblrow

    OUT0("Starting TpetraCrsGraph_insertGlobalIndices!")

    ! Duplicate Tpetra_CrsGraph_CreateTestGraph_A here because
    ! we want insertGlobalIndices to appear explicitly.
    rmap = TpetraMap(Tpetra_GLOBAL_INVALID, test_graph_num_row(), comm)
    Graph = TpetraCrsGraph(rmap, test_graph_max_entries_per_row(), TpetraStaticProfile)

    do irow=1,test_graph_num_row()
      call Tpetra_CrsGraph_GetTestGraphRow_A(comm, irow, gblrow, cols_expected, nnz)
      call Graph%insertGlobalIndices(gblrow, cols_expected)
      deallocate(cols_expected)
    enddo

    !---------------------------------------------------------------
    ! CHECK results
    !---------------------------------------------------------------
    do irow=1,test_graph_num_row()
      call Tpetra_CrsGraph_GetTestGraphRow_A(comm, irow, gblrow, cols_expected, nnz)
      allocate(cols(nnz))
      call Graph%getGlobalRowCopy(gblrow, cols, nument)
      TEST_EQUALITY(nnz, int(nument))
      TEST_ARRAY_EQUALITY(cols, cols_expected)
      deallocate(cols_expected)
      deallocate(cols)
    enddo

    call Graph%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_insertGlobalIndices!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_insertGlobalIndices)

  ! ----------------------------removeLocalIndices---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_removeLocalIndices)
    type(TpetraCrsGraph) :: Graph
    type(TpetraMap) :: rowMap, colMap
    integer(int_type) :: lclNumRows, lclNumCols, lclRow
    integer(global_ordinal_type) gblNumRows
    integer(int_type) :: numProcs
    integer(int_type), dimension(:), allocatable :: lclColInds
    integer(int_type), dimension(:), allocatable :: curLclColInds
    integer(size_type), parameter :: max_entries_per_row=3
    integer(size_type) :: numEnt=0

    OUT0("Starting TpetraCrsGraph_removeLocalIndices!")

    numProcs = comm%getSize()
    lclNumRows = 10
    gblNumRows = lclNumRows * numProcs
    rowMap = TpetraMap(gblNumRows, lclNumRows, comm)
    lclNumCols = lclNumRows
    colMap = rowMap
    Graph = TpetraCrsGraph(rowMap, colMap, max_entries_per_row, TpetraStaticProfile)

    allocate(lclColInds(2))

    do lclRow = 1, lclNumRows
      lclColInds(1) = 1 + modulo(lclRow + 0, lclNumCols)
      lclColInds(2) = 1 + modulo(lclRow + 1, lclNumCols)
      TEST_NOTHROW(call Graph%insertLocalIndices(lclRow, lclColInds))
    end do

    ! Now remove just row 5
    lclRow = 5
    TEST_NOTHROW(call Graph%removeLocalIndices(lclRow))

    ! Make the arrays bigger than necessary, just to make sure that
    ! the methods behave correctly.
    allocate(curLclColInds(5))

    do lclRow = 1, lclNumRows
      TEST_NOTHROW(call Graph%getLocalRowCopy(lclRow, curLclColInds, numEnt))

      ! We should have this row only to be empty
      if(lclRow .EQ. 5) then
        TEST_EQUALITY(int(numEnt), 0)
      else
        TEST_EQUALITY(int(numEnt), 2)

      TEST_EQUALITY(curLclColInds(1), 1 + modulo(lclRow + 0, lclNumCols));
      TEST_EQUALITY(curLclColInds(2), 1 + modulo(lclRow + 1, lclNumCols));
      endif

    enddo

    call Graph%release(); TEST_IERR()
    call rowMap%release(); TEST_IERR()
    call colMap%release(); TEST_IERR()
    deallocate(lclColInds)
    deallocate(curLclColInds)

    OUT0("Finished TpetraCrsGraph_removeLocalIndices!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_removeLocalIndices)

  ! ---------------------------------getComm---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getComm)
    type(TpetraCrsGraph) :: Graph
    type(TeuchosComm) :: gcomm

    OUT0("Starting TpetraCrsGraph_getComm!")

    call Tpetra_CrsGraph_CreateBasic(comm, Graph)

    gcomm = Graph%getComm(); TEST_IERR()

    ! We only test the comm returned has the same rank and size.  More
    ! comprehensive testing is (assumed to be) done in Teuchos itself.
    TEST_ASSERT(gcomm%getRank() == comm%getRank())
    TEST_ASSERT(gcomm%getSize() == comm%getSize())

    call gcomm%release(); TEST_IERR()
    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getComm!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getComm)

  ! --------------------------------getRowMap--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getRowMap)
    type(TpetraMap) :: rmap, tmpmap
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_getRowMap!")

    ! Consistent with Tpetra_CrsGraph_CreateTestGraph_B
    rmap = TpetraMap(TPETRA_GLOBAL_INVALID, test_graph_num_local(), comm)

    call Tpetra_CrsGraph_CreateTestGraph_B(comm,Graph)
    call Graph%fillComplete()

    tmpmap = Graph%getRowMap()

    TEST_EQUALITY(tmpmap%getNodeNumElements(), rmap%getNodeNumElements())
    TEST_ASSERT(tmpmap%isSameAs(rmap))

    call tmpmap%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()
    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getRowMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getRowMap)

  ! --------------------------------getColMap--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getColMap)
    type(TpetraMap) :: cmap, tmpmap
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_getColMap!")

    ! Consistent with Tpetra_CrsGraph_CreateTestGraph_B
    cmap = TpetraMap(TPETRA_GLOBAL_INVALID, test_graph_num_local(), comm)

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()
    tmpmap = Graph%getColMap()

    TEST_EQUALITY(tmpmap%getNodeNumElements(), cmap%getNodeNumElements())
    TEST_ASSERT(tmpmap%isSameAs(cmap))

    call tmpmap%release(); TEST_IERR()
    call Graph%release(); TEST_IERR();
    call cmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getColMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getColMap)

  ! -------------------------------getDomainMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getDomainMap)
    type(TpetraMap) :: Map, tmpmap
    type(TpetraCrsGraph) :: Graph
    integer(global_ordinal_type) :: num_global

    OUT0("Starting TpetraCrsGraph_getDomainMap!")

    ! create Map and Graph
    num_global = test_graph_num_local()*comm%getSize()
    ! Consistent with Tpetra_CrsGraph_CreateTestGraph_B
    Map = TpetraMap(num_global, comm); TEST_IERR()
    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    tmpmap = Graph%getDomainMap(); TEST_IERR()
    TEST_ASSERT(tmpmap%isSameAs(Map))

    call tmpmap%release(); TEST_IERR()
    call Graph%release(); TEST_IERR()
    call Map%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getDomainMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getDomainMap)

  ! -------------------------------getRangeMap-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getRangeMap)
    type(TpetraMap) :: Map, tmpmap
    type(TpetraCrsGraph) :: Graph
    integer(size_type), parameter :: izero=0, ione=1
    integer(global_ordinal_type) :: num_global

    OUT0("Starting TpetraCrsGraph_getRangeMap!")

    ! create Map and Graph
    num_global = test_graph_num_local()*comm%getSize()
    ! Consistent with Tpetra_CrsGraph_CreateTestGraph_B
    Map = TpetraMap(num_global, comm); TEST_IERR()
    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    tmpmap = Graph%getRangeMap(); TEST_IERR()
    TEST_ASSERT(tmpmap%isSameAs(Map))

    call tmpmap%release(); TEST_IERR()
    call Graph%release(); TEST_IERR()
    call Map%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getRangeMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getRangeMap)

  ! -------------------------------getImporter-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getImporter)
    type(TpetraMap) :: rmap, cmap, imprmap, impcmap
    type(TpetraCrsGraph) :: Graph
    type(TpetraImport) :: importer
    integer(size_type) :: num_ent_per_row

    OUT0("Starting TpetraCrsGraph_getImporter!")

    rmap = TpetraMap(TPETRA_GLOBAL_INVALID, test_graph_num_local(), comm); TEST_IERR()
    cmap = TpetraMap(TPETRA_GLOBAL_INVALID, 2*test_graph_num_local(), comm); TEST_IERR()
    ! must allocate enough for all submitted indices.
    num_ent_per_row = 2
    Graph = TpetraCrsGraph(rmap, cmap, num_ent_per_row)
    call Graph%fillComplete()

    importer = Graph%getImporter()
    imprmap = importer%getSourceMap()
    TEST_ASSERT(rmap%isSameAs(imprmap))
    impcmap = importer%getTargetMap()
    TEST_ASSERT(cmap%isSameAs(impcmap))

    call imprmap%release(); TEST_IERR()
    call impcmap%release(); TEST_IERR()
    call importer%release(); TEST_IERR()
    call Graph%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()
    call cmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getImporter!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getImporter)

  ! -------------------------------getExporter-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getExporter)
    type(TpetraMap) :: rmap, cmap, expmap
    type(TpetraCrsGraph) :: Graph
    type(TpetraExport) :: exporter
    integer(size_type) :: num_ent_per_row

    OUT0("Starting TpetraCrsGraph_getExporter!")

    rmap = TpetraMap(TPETRA_GLOBAL_INVALID, test_graph_num_local(), comm); TEST_IERR()
    cmap = TpetraMap(TPETRA_GLOBAL_INVALID, 2*test_graph_num_local(), comm); TEST_IERR()
    ! must allocate enough for all submitted indices.
    num_ent_per_row = 2
    Graph = TpetraCrsGraph(rmap, cmap, num_ent_per_row)

    ! make range map different so exporter will exist
    call Graph%fillComplete(rmap, cmap)

    exporter = Graph%getExporter()
    expmap = exporter%getSourceMap()
    TEST_ASSERT(rmap%isSameAs(expmap))

    call exporter%release(); TEST_IERR()
    call Graph%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()
    call cmap%release(); TEST_IERR()
    call expmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getExporter!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getExporter)

  ! -----------------------------getGlobalNumRows----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalNumRows)
    type(TpetraMap) :: rmap, cmap
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: num_ent_per_row, ires
    integer :: num_local
    integer :: commsize

    OUT0("Starting TpetraCrsGraph_getGlobalNumRows!")

    num_local = 1
    rmap = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    cmap = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    ! must allocate enough for all submitted indices.
    num_ent_per_row = 2
    Graph = TpetraCrsGraph(rmap, cmap, num_ent_per_row, TpetraStaticProfile)

    commsize = comm%getSize()
    ires = Graph%getGlobalNumRows()
    TEST_EQUALITY(ires, int(commsize,size_type))

    call Graph%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()
    call cmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getGlobalNumRows!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalNumRows)

  ! -----------------------------getGlobalNumCols----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalNumCols)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsGraph_getGlobalNumCols!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    ires = Graph%getGlobalNumCols()
    TEST_EQUALITY(ires, int(test_graph_num_local()*comm%getSize(),size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getGlobalNumCols!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalNumCols)

  ! ---------------------------getGlobalNumEntries---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalNumEntries)
    type(TpetraMap) :: rmap, cmap
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsGraph_getGlobalNumEntries!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    ires = Graph%getGlobalNumEntries()
    ! match num_ent_per_row in CreateTestGraph
    TEST_EQUALITY(ires, int(2*comm%getSize(),size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getGlobalNumEntries!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalNumEntries)

  ! ------------------------getGlobalMaxNumRowEntries------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalMaxNumRowEntries)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) ::  ires

    OUT0("Starting TpetraCrsGraph_getGlobalMaxNumRowEntries!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    ires = Graph%getGlobalMaxNumRowEntries()
    ! match num_ent_per_row in CreateTestGraph
    TEST_EQUALITY(ires, int(2,size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getGlobalMaxNumRowEntries!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalMaxNumRowEntries)

  ! ------------------------------getNodeNumRows------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeNumRows)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) ::  ires

    OUT0("Starting TpetraCrsGraph_getNodeNumRows!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    ires = Graph%getNodeNumRows()
    TEST_EQUALITY(ires, int(test_graph_num_local(), size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNodeNumRows!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeNumRows)

  ! ------------------------------getNodeNumCols------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeNumCols)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsGraph_getNodeNumCols!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    ires = Graph%getNodeNumCols()
    TEST_EQUALITY(ires, int(test_graph_num_local(), size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNodeNumCols!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeNumCols)

  ! ----------------------------getNodeNumEntries----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeNumEntries)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsGraph_getNodeNumEntries!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    ires = Graph%getNodeNumEntries()
    TEST_EQUALITY(ires, int(2,size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNodeNumEntries!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeNumEntries)

  ! -------------------------getNodeMaxNumRowEntries-------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeMaxNumRowEntries)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsGraph_getNodeMaxNumRowEntries!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    ires = Graph%getNodeMaxNumRowEntries()
    TEST_EQUALITY(ires, int(2,size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNodeMaxNumRowEntries!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeMaxNumRowEntries)

  ! -------------------------getNumEntriesInGlobalRow------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNumEntriesInGlobalRow)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires
    type(TpetraMap) :: rmap
    integer(global_ordinal_type) :: myrowind
    integer :: lclrow

    OUT0("Starting TpetraCrsGraph_getNumEntriesInGlobalRow!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    lclrow = 1
    rmap = Graph%getRowMap()
    myrowind = rmap%getGlobalElement(lclrow);
    ires = Graph%getNumEntriesInGlobalRow(myrowind)
    TEST_EQUALITY(ires, int(2,size_type))

    call Graph%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNumEntriesInGlobalRow!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNumEntriesInGlobalRow)

  ! -------------------------getNumEntriesInLocalRow-------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNumEntriesInLocalRow)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires
    integer :: lclrow

    OUT0("Starting TpetraCrsGraph_getNumEntriesInLocalRow!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    lclrow = 1
    ires = Graph%getNumEntriesInLocalRow(lclrow)
    TEST_EQUALITY(ires, int(2,size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNumEntriesInLocalRow!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNumEntriesInLocalRow)

  ! --------------------------getNodeAllocationSize--------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeAllocationSize)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsGraph_getNodeAllocationSize!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    ires = Graph%getNodeAllocationSize()
    ! match num_ent_per_row in CreateTestGraph
    TEST_EQUALITY(ires, int(2,size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNodeAllocationSize!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeAllocationSize)

  ! --------------------getNumAllocatedEntriesInGlobalRow--------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNumAllocatedEntriesInGlobalRow)
    type(TpetraMap) :: rmap
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires
    integer(global_ordinal_type) :: myrowind
    integer :: lclrow

    OUT0("Starting TpetraCrsGraph_getNumAllocatedEntriesInGlobalRow!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    lclrow = 1
    rmap = Graph%getRowMap()
    myrowind = rmap%getGlobalElement(lclrow);
    ires = Graph%getNumAllocatedEntriesInGlobalRow(myrowind)
    TEST_EQUALITY(ires, int(2,size_type))

    call Graph%release(); TEST_IERR()
    call rmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNumAllocatedEntriesInGlobalRow!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNumAllocatedEntriesInGlobalRow)

  ! ---------------------getNumAllocatedEntriesInLocalRow--------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNumAllocatedEntriesInLocalRow)
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: ires
    integer :: lclrow

    OUT0("Starting TpetraCrsGraph_getNumAllocatedEntriesInLocalRow!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    lclrow = 1
    ires = Graph%getNumAllocatedEntriesInLocalRow(lclrow)
    TEST_EQUALITY(ires, int(2,size_type))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNumAllocatedEntriesInLocalRow!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNumAllocatedEntriesInLocalRow)

  ! ------------------------------getProfileType------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getProfileType)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_getProfileType!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    TEST_ASSERT(graph%getProfileType() == TpetraStaticProfile)

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getProfileType!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getProfileType)

  ! -----------------------------getGlobalRowCopy----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalRowCopy)
    type(TpetraCrsGraph) :: Graph
    integer :: irow, nnz
    integer(size_type) :: nument;
    integer(global_ordinal_type), allocatable :: cols_expected(:)
    integer(global_ordinal_type), allocatable :: cols(:)
    integer(global_ordinal_type) :: gblrow

    OUT0("Starting TpetraCrsGraph_getGlobalRowCopy!")

    call Tpetra_CrsGraph_CreateTestGraph_A(comm, Graph)
    do irow=1,test_graph_num_row()
      call Tpetra_CrsGraph_GetTestGraphRow_A(comm, irow, gblrow, cols_expected, nnz)
      allocate(cols(nnz))
      call Graph%getGlobalRowCopy(gblrow, cols, nument)
      TEST_EQUALITY(nnz, int(nument))
      TEST_ARRAY_EQUALITY(cols, cols_expected)
      deallocate(cols_expected)
      deallocate(cols)
    enddo

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getGlobalRowCopy!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getGlobalRowCopy)

  ! -----------------------------getLocalRowCopy------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getLocalRowCopy)
    type(TpetraCrsGraph) :: Graph
    type(TpetraMap) :: rowMap, colMap
    integer(int_type) :: lclNumRows, lclNumCols, lclRow
    integer(global_ordinal_type) gblNumRows
    integer(int_type) :: numProcs
    integer(int_type), dimension(:), allocatable :: lclColInds
    integer(int_type), dimension(:), allocatable :: curLclColInds
    integer(size_type), parameter :: max_entries_per_row=3
    integer(size_type) :: numEnt=0

    OUT0("Finished TpetraCrsGraph_insertLocalIndices!")

    numProcs = comm%getSize()
    lclNumRows = 10
    gblNumRows = lclNumRows * numProcs
    rowMap = TpetraMap(gblNumRows, lclNumRows, comm)
    lclNumCols = lclNumRows
    colMap = rowMap
    Graph = TpetraCrsGraph(rowMap, colMap, max_entries_per_row, TpetraStaticProfile)

    allocate(lclColInds(2))

    do lclRow = 1, lclNumRows
      lclColInds(1) = 1 + modulo(lclRow + 0, lclNumCols)
      lclColInds(2) = 1 + modulo(lclRow + 1, lclNumCols)
      TEST_NOTHROW(call Graph%insertLocalIndices(lclRow, lclColInds))
    end do

    ! Make the arrays bigger than necessary, just to make sure that
    ! the methods behave correctly.
    allocate(curLclColInds(5))

    do lclRow = 1, lclNumRows
      TEST_NOTHROW(call Graph%getLocalRowCopy(lclRow, curLclColInds, numEnt))
      TEST_EQUALITY(int(numEnt), 2)

      TEST_EQUALITY(curLclColInds(1), 1 + modulo(lclRow + 0, lclNumCols));
      TEST_EQUALITY(curLclColInds(2), 1 + modulo(lclRow + 1, lclNumCols));
    enddo

    call Graph%release(); TEST_IERR()
    call rowMap%release(); TEST_IERR()
    call colMap%release(); TEST_IERR()
    deallocate(lclColInds)
    deallocate(curLclColInds)

    OUT0("Finished TpetraCrsGraph_getLocalRowCopy!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getLocalRowCopy)

  ! ------------------------------getNodeRowPtrs------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeRowPtrs)
    type(TpetraCrsGraph) :: Graph
    integer(size_type), allocatable :: rowpointers(:)

    OUT0("Starting TpetraCrsGraph_getNodeRowPtrs!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    allocate(rowpointers(1))
    TEST_THROW(call Graph%getNodeRowPtrs(rowpointers(:)))

    deallocate(rowpointers)
    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNodeRowPtrs!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodeRowPtrs)

  ! ---------------------------getNodePackedIndices--------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodePackedIndices)
    type(TpetraCrsGraph) :: Graph
    integer(size_type), allocatable :: columnindices(:)

    OUT0("Starting TpetraCrsGraph_getNodePackedIndices!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    allocate(columnindices(1))
    TEST_THROW(call Graph%getNodePackedIndices(columnindices(:)))

    deallocate(columnindices)
    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_getNodePackedIndices!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_getNodePackedIndices)

  ! --------------------------------hasColMap--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_hasColMap)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_hasColMap!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    TEST_ASSERT(Graph%hasColMap())

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_hasColMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_hasColMap)

  ! -----------------------------isLocallyIndexed----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isLocallyIndexed)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_isLocallyIndexed!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()

    TEST_ASSERT(Graph%isLocallyIndexed())

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_isLocallyIndexed!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isLocallyIndexed)

  ! ----------------------------isGloballyIndexed----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isGloballyIndexed)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_isGloballyIndexed!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    TEST_ASSERT(Graph%isGloballyIndexed())

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_isGloballyIndexed!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isGloballyIndexed)

  ! ------------------------------isFillComplete------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isFillComplete)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_isFillComplete!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    TEST_ASSERT(.not. Graph%isFillComplete())
    call Graph%fillComplete()
    TEST_ASSERT(Graph%isFillComplete())

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_isFillComplete!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isFillComplete)

  ! -------------------------------isFillActive------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isFillActive)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_isFillActive!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    TEST_ASSERT(Graph%isFillActive())
    call Graph%fillComplete()
    TEST_ASSERT(.not. Graph%isFillActive())

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_isFillActive!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isFillActive)

  ! ---------------------------------isSorted--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isSorted)
    type(TpetraMap) :: map
    type(TpetraCrsGraph) :: Graph
    integer(size_type) :: nument
    integer :: num_local

    OUT0("Starting TpetraCrsGraph_isSorted!")

    ! create a Map
    num_local = 10;
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm)

    nument = 4
    Graph = TpetraCrsGraph(map, map, nument); TEST_IERR()
    TEST_ASSERT(graph%isSorted())

    call Graph%release(); TEST_IERR()
    call map%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_isSorted!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isSorted)

  ! ----------------------------isStorageOptimized---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isStorageOptimized)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_isStorageOptimized!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    TEST_ASSERT(.not. Graph%isStorageOptimized())
    call Graph%fillComplete()
    TEST_ASSERT(Graph%isStorageOptimized())

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_isStorageOptimized!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_isStorageOptimized)

  ! -----------------------------supportsRowViews----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_supportsRowViews)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_supportsRowViews!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    TEST_ASSERT(Graph%supportsRowViews())

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_supportsRowViews!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_supportsRowViews)

  ! -------------------------------description-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_description)
    type(TpetraCrsGraph) :: Graph
    character(kind=C_CHAR) :: fresult

    OUT0("Starting TpetraCrsGraph_description!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    TEST_NOTHROW(fresult = Graph%description())

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_description!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_description)

  ! ------------------------------replaceColMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_replaceColMap)
    type(TpetraMap) :: Map, newcolmap
    type(TpetraExport) :: fresult
    type(TpetraCrsGraph) :: Graph
    integer(global_ordinal_type) :: num_global

    OUT0("Starting TpetraCrsGraph_replaceColMap!")

    num_global = 4*comm%getSize()
    Map = TpetraMap(num_global, comm); TEST_IERR()
    Graph = TpetraCrsGraph(Map, Map, 1_size_type, TpetraStaticProfile)

    ! Make a new colmap
    newcolmap = TpetraMap(num_global, comm); TEST_IERR()

    TEST_NOTHROW(call Graph%replaceColMap(newcolmap))

    call Graph%release(); TEST_IERR()
    call Map%release(); TEST_IERR()
    call newcolmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_replaceColMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_replaceColMap)

  ! -----------------------replaceDomainMapAndImporter------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_replaceDomainMapAndImporter)
    type(TpetraMap) :: Map, newdomainmap, newcolmap
    type(TpetraImport) :: newimporter
    type(TpetraExport) :: fresult
    type(TpetraCrsGraph) :: Graph
    integer(size_type), parameter :: izero=0, ione=1
    integer(global_ordinal_type) :: num_global

    OUT0("Starting TpetraCrsGraph_replaceDomainMapAndImporter!")

    num_global = 4*comm%getSize()
    Map = TpetraMap(num_global, comm); TEST_IERR()
    Graph = TpetraCrsGraph(Map, Map, izero, TpetraStaticProfile)
    call Graph%fillComplete()     ! Needed for gloStaticstants

    ! Make a new colmap
    num_global = 5*comm%getSize()
    newdomainmap = TpetraMap(num_global, comm); TEST_IERR()
    newcolmap =  Graph%getColMap(); TEST_IERR()
    newimporter = TpetraImport(newdomainmap, newcolmap ); TEST_IERR()

    TEST_NOTHROW(call Graph%replaceDomainMapAndImporter(newdomainmap, newimporter))

    call newimporter%release(); TEST_IERR()
    call newdomainmap%release(); TEST_IERR()
    call newcolmap%release(); TEST_IERR()
    call Graph%release(); TEST_IERR()
    call Map%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_replaceDomainMapAndImporter!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_replaceDomainMapAndImporter)

  ! -----------------------removeEmptyProcessesInPlace------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_removeEmptyProcessesInPlace)
    type(TpetraMap) :: Map
    type(TpetraCrsGraph) :: Graph
    integer(size_type), parameter :: izero=0, ione=1
    integer(global_ordinal_type) :: num_global

    OUT0("Starting TpetraCrsGraph_removeEmptyProcessesInPlace!")

    ! create Map and Graph
    num_global = 4*comm%getSize()
    Map = TpetraMap(num_global, comm); TEST_IERR()
    Graph = TpetraCrsGraph(Map, Map, izero, TpetraStaticProfile)

    TEST_NOTHROW(call Graph%removeEmptyProcessesInPlace(Map))

    call Graph%release(); TEST_IERR()
    call Map%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_removeEmptyProcessesInPlace!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_removeEmptyProcessesInPlace)

  ! ---------------------------haveGlobalConstants---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_haveGlobalConstants)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_haveGlobalConstants!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)
    call Graph%fillComplete()     ! Needed for global constants

    TEST_ASSERT(Graph%haveGlobalConstants())

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_haveGlobalConstants!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_haveGlobalConstants)

  ! --------------------------computeGlobalConstants-------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsGraph_computeGlobalConstants)
    type(TpetraCrsGraph) :: Graph

    OUT0("Starting TpetraCrsGraph_computeGlobalConstants!")

    call Tpetra_CrsGraph_CreateTestGraph_B(comm, Graph)

    TEST_NOTHROW(call Graph%computeGlobalConstants(.true.))

    call Graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsGraph_computeGlobalConstants!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsGraph_computeGlobalConstants)

end program test_tpetraCrsGraph
