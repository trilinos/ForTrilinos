! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
!
! Routines for creating matrices for testing
! Functions would be more elegant, but rcp errors occurred
!
#include "ForTrilinos_config.h"
module test_Tpetra_crsgraph_helper
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use forTpetra

  implicit none

contains
  ! helper functions to define simple tests
  function test_graph_num_local() result (n)
    integer :: n
    n = 10
  end function test_graph_num_local

  function test_graph_num_row() result (n)
    integer :: n
    n = 4
  end function test_graph_num_row

  function test_graph_max_entries_per_row() result (n)
    integer(size_type) :: n
    n = 3
  end function test_graph_max_entries_per_row

  ! ----------------------- Create Basic Graph ------------------------------- !
  !   do not call fillComplete
  subroutine Tpetra_CrsGraph_CreateBasic(comm,Graph)
    type(TeuchosComm), intent(in) :: comm
    type(TpetraCrsGraph), intent(out) :: Graph
    type(TpetraMap) :: map

    ! create a Map
    Map = TpetraMap(TPETRA_GLOBAL_INVALID, test_graph_num_local(), comm)
    Graph = TpetraCrsGraph(Map, Map, int(test_graph_num_local(),size_type), TpetraStaticProfile)

    call map%release()
    return
  end subroutine Tpetra_CrsGraph_CreateBasic
  ! ----------------------- Create Test Graph ------------------------------- !
  subroutine Tpetra_CrsGraph_CreateTestGraph_B(comm,Graph)
    type(TeuchosComm), intent(in) :: comm
    type(TpetraCrsGraph), intent(out) :: Graph
    type(TpetraMap) :: rmap, cmap
    integer(size_type) :: num_ent_per_row
    integer(global_ordinal_type) :: myrowind
    integer(global_ordinal_type), allocatable :: indices(:)
    integer :: lclrow

    ! create maps
    rmap = TpetraMap(TPETRA_GLOBAL_INVALID, test_graph_num_local(), comm)
    cmap = TpetraMap(TPETRA_GLOBAL_INVALID, test_graph_num_local(), comm)

    ! must allocate enough for all submitted indices.
    num_ent_per_row = 2
    Graph = TpetraCrsGraph(rmap, cmap, num_ent_per_row, TpetraStaticProfile)

    lclrow = 1
    myrowind = rmap%getGlobalElement(lclrow);

    ! insertGlobalIndices doesn't do column Map filtering, so we have to test
    ! whether each of the column indices to insert is invalid.
    if (cmap%isNodeGlobalElement(myrowind)) then
      if (cmap%isNodeGlobalElement(myrowind+1)) then
        allocate(indices(2))
        indices = [myrowind, myrowind+1]
        call Graph%insertGlobalIndices(myrowind, indices)
        deallocate(indices)
      else
        allocate(indices(1))
        indices(1) = myrowind
        call Graph%insertGlobalIndices(myrowind, indices)
        deallocate(indices)
      end if
    end if

    call rmap%release();
    call cmap%release();

    return
  end subroutine Tpetra_CrsGraph_CreateTestGraph_B
  ! ------------------------------- Create Test Graph2 ------------------------- !
  !   do not call fillComplete
  !  This is the graph of the following matrix that is created in
  !     test_tpetra_crsmatrix_helper
  !  [2 1           ]
  !  [1 1 1         ]
  !  [  1 1 1       ]
  !  [   . . .      ]
  !  [     . . .    ]
  !  [       . . .  ]
  !  [         1 1 1]
  !  [           1 2]
  ! this matrix has an eigenvalue lambda=3, with eigenvector v = [1 ... 1]
  subroutine Tpetra_CrsGraph_GetTestGraphRow_A(comm, irow, gblrow, cols, nnz)
    type(TeuchosComm), intent(in) :: comm
    integer, intent(in) :: irow
    integer(global_ordinal_type), intent(out), allocatable :: cols(:)
    integer, intent(out) :: nnz
    integer(global_ordinal_type), intent(out) :: gblrow
    integer :: num_images, my_image_id
    num_images = comm%getSize()
    my_image_id = comm%getRank()

    if (num_images < 2) then
      gblrow = irow
    else
      gblrow = my_image_id*test_graph_num_row() + irow
    endif
    if (gblrow == 1) then
      nnz = 2
      allocate(cols(nnz))
      cols = [gblrow, gblrow+1]
    else if (gblrow == num_images*test_graph_num_row()) then
      nnz = 2
      allocate(cols(nnz))
      cols = [gblrow-1, gblrow]
    else
      nnz = 3
      allocate(cols(nnz))
      cols = [gblrow-1, gblrow, gblrow+1]
    end if

  end subroutine Tpetra_CrsGraph_GetTestGraphRow_A

  subroutine Tpetra_CrsGraph_CreateTestGraph_A(comm,Graph)

    type(TeuchosComm), intent(in) :: comm
    type(TpetraCrsGraph), intent(out) :: Graph
    type(TpetraMap) :: rmap
    integer :: irow, nnz
    integer(global_ordinal_type), allocatable :: cols(:)
    integer(global_ordinal_type) :: gblrow

    rmap = TpetraMap(Tpetra_GLOBAL_INVALID, test_graph_num_row(), comm)
    Graph = TpetraCrsGraph(rmap, test_graph_max_entries_per_row(), TpetraStaticProfile)

    do irow=1,test_graph_num_row()
      call Tpetra_CrsGraph_GetTestGraphRow_A(comm, irow, gblrow, cols, nnz)
      call Graph%insertGlobalIndices(gblrow, cols)
      deallocate(cols)
    enddo

    call rmap%release();

    return
  end subroutine Tpetra_CrsGraph_CreateTestGraph_A
  ! ------------------------------- Get Test Graph2 Results ------------------ !
  ! Hard-code the results graph for getAllValues test from the TestGraph2
  ! routine above

! MDM - commented this out because we don't currently use it
! TO add back in when it's actually tested or delete this if not needed.

!  subroutine Tpetra_CrsGraph_getTestGraphResults(comm, row_ptrs, cols)
!    type(TeuchosComm), intent(in) :: comm
!    integer(size_type), dimension(:), intent(out) :: row_ptrs
!    integer(int_type), dimension(:), intent(out) :: cols
!    integer :: nnz, irow, idx, i
!    integer(global_ordinal_type), allocatable :: cols_row(:)
!    integer(global_ordinal_type) :: gblrow
!
!    row_ptrs(1)=1
!    do irow=1,test_graph_num_row()
!      idx=row_ptrs(irow)
!      call Tpetra_CrsGraph_GetTestGraphRow(comm, irow, gblrow, cols_row, nnz)
!      do i=1,nnz
!        cols(idx+i-1) = cols_row(i)
!      enddo
!      row_ptrs(irow+1)=row_ptrs(irow)+nnz
!      deallocate(cols_row)
!    enddo
!
!    return
!  end subroutine Tpetra_CrsGraph_getTestGraphResults
end module test_tpetra_crsgraph_helper
