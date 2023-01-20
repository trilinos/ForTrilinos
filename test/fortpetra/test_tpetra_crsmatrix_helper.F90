! Copyright 2017-2018, UT-Batte+nnzlle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
!
! Routines for creating matrices for testing
! Functions would be more elegant, but rcp errors occurred
!
#include "ForTrilinos_config.h"
module test_Tpetra_crsmatrix_helper
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  implicit none

contains
  ! helper functions to define simple tests
  function test_matrix_num_local() result (n)
    integer :: n
    n = 10
  end function test_matrix_num_local

  function test_matrix_num_row() result (n)
    integer :: n
    n = 4
  end function test_matrix_num_row

  function test_matrix_max_entries_per_row() result (n)
    integer(size_type) :: n
    n = 3
  end function test_matrix_max_entries_per_row

  ! ----------------------- Create Identity Matrix ----------------------------- !
  !   do not call fillComplete
  subroutine Tpetra_CrsMatrix_CreateIdentity(comm,Mat)
    type(TeuchosComm), intent(in) :: comm
    type(TpetraCrsMatrix), intent(out) :: Mat
    type(TpetraMap) :: map
    integer :: irow
    integer(global_ordinal_type) :: base, gblrow, cols(1)
    real(mag_type) :: vals(1)=[1.]

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, test_matrix_num_local(), comm)

    ! create the identity matrix
    base = test_matrix_num_local() * comm%getRank()
    Mat = TpetraCrsMatrix(map, map, 1_size_type)
    do irow = 1, test_matrix_num_local()
      gblrow = base + int(irow, kind=global_ordinal_type)
      cols(1) = gblrow
      vals(1) = 1.d0
      call Mat%insertGlobalValues(gblrow, cols, vals)
    end do

    call map%release()
    return
  end subroutine Tpetra_CrsMatrix_CreateIdentity

  ! ------------------------------- Create Test Matrix ------------------------- !
  !   do not call fillComplete
  !  create the following matrix:
  !  [2 1           ]
  !  [1 1 1         ]
  !  [  1 1 1       ]
  !  [   . . .      ]
  !  [     . . .    ]
  !  [       . . .  ]
  !  [         1 1 1]
  !  [           1 2]
  ! this matrix has an eigenvalue lambda=3, with eigenvector v = [1 ... 1]
  subroutine Tpetra_CrsMatrix_CreateTestMatrix_A(comm,Mat)
    type(TeuchosComm), intent(in) :: comm
    type(TpetraCrsMatrix), intent(out) :: Mat
    type(TpetraMap) :: map
    integer(global_ordinal_type), allocatable :: cols(:)
    real(scalar_type), allocatable :: vals(:)
    integer :: irow, nnz
    integer(global_ordinal_type) :: gblrow

    ! create a Map and matrix
    map = TpetraMap(TPETRA_GLOBAL_INVALID, test_matrix_num_row(), comm)

    ! NOTE -  could consider doing a 2D distribution to teach that
    !rowmap = TpetraMap(TPETRA_GLOBAL_INVALID, test_matrix_num_row(), comm)
    !nnz=test_matrix_num_row()()*test_matrix_max_entries_per_row()
    !colmap = TpetraMap(TPETRA_GLOBAL_INVALID, nnz, comm)
    !Mat = TpetraCrsMatrix(rowmap, colmap, test_matrix_max_entries_per_row())
    Mat = TpetraCrsMatrix(map, test_matrix_max_entries_per_row())

    do irow=1,test_matrix_num_row()
      call Tpetra_CrsMatrix_GetTestMatrixRow_A(comm, irow, gblrow, cols, vals, nnz)
      call Mat%insertGlobalValues(gblrow, cols, vals)
      deallocate(cols,vals)
    enddo

    call map%release()

  return
  end subroutine Tpetra_CrsMatrix_CreateTestMatrix_A

  ! ------------------------------- Get Test Matrix Results ------------------ !
  ! Hard-code the results matrices for getAllValues test from the TestMatrix
  ! routine above
  subroutine Tpetra_CrsMatrix_getTestMatrixResults(comm, row_ptrs, cols, vals)
    type(TeuchosComm), intent(in) :: comm
    integer(size_type), dimension(:), intent(out) :: row_ptrs
    integer(int_type), dimension(:), intent(out) :: cols
    real(scalar_type), dimension(:), intent(out) :: vals
    integer(global_ordinal_type), dimension(:), allocatable :: cols_row
    real(scalar_type), dimension(:), allocatable :: vals_row
    integer :: num_images, my_image_id, nnz, irow, idx, i
    integer(global_ordinal_type) :: gblrow

    num_images = comm%getSize()
    my_image_id = comm%getRank()

    ! insertMatValues routine above uses global indices for row and cols.
    ! here we return local row but global indices for column
    row_ptrs(1)=1
    do irow=1,test_matrix_num_row()
      idx=row_ptrs(irow)
      call Tpetra_CrsMatrix_GetTestMatrixRow_A(comm, irow, gblrow, cols_row, vals_row, nnz)
      do i=1,nnz
        cols(idx+i-1) = cols_row(i)
        vals(idx+i-1) = vals_row(i)
      enddo
      deallocate(cols_row,vals_row)
      row_ptrs(irow+1)=row_ptrs(irow)+nnz
    enddo

    return
  end subroutine Tpetra_CrsMatrix_getTestMatrixResults
  ! ------------------------------- GetTestMatrixRow   ------------------------- !
  ! Utility subroutine to handle the allocation logic for  creating the testMatrix
  !   and the expect results:
  ! Tpetra_CrsMatrix_CreateTestMatrix and  Tpetra_CrsMatrix_getTestMatrixResults
  !
  subroutine Tpetra_CrsMatrix_GetTestMatrixRow_A(comm, irow, gblrow, cols, vals, nnz)
    type(TeuchosComm), intent(in) :: comm
    integer, intent(in) :: irow
    integer(global_ordinal_type), intent(out) :: gblrow
    integer(global_ordinal_type), intent(out), allocatable :: cols(:)
    real(scalar_type), intent(out), allocatable :: vals(:)
    integer, intent(out) :: nnz
    integer :: num_images, my_image_id

    num_images = comm%getSize()
    my_image_id = comm%getRank()

    if (num_images < 2) then
      gblrow = irow
    else
      gblrow = my_image_id*test_matrix_num_row() + irow
    endif
    if (gblrow == 1) then
      nnz = 2
      allocate(cols(nnz)); allocate(vals(nnz));
      cols(1:nnz) = [gblrow, gblrow+1]
      vals(1:nnz) = [2.d0, 1.d0]
    else if (gblrow == num_images*test_matrix_num_row()) then
      nnz = 2;
      allocate(cols(nnz)); allocate(vals(nnz));
      cols(1:nnz) = [gblrow-1, gblrow]
      vals(1:nnz) = [1.d0, 2.d0]
    else
      nnz = 3;
      allocate(cols(nnz)); allocate(vals(nnz));
      cols = [gblrow-1, gblrow, gblrow+1]
      vals = [1.d0, 1.d0, 1.d0]
    end if

    return
  end subroutine Tpetra_CrsMatrix_GetTestMatrixRow_A
  ! ------------------------------- Create Labelled Test Matrix ------------------------- !
  !  do not call fillComplete
  !  Create a matrix where things are globally labelled in a way to allow easy debugging
  !  [2 3           ]
  !  [  3 4         ]
  !  [    4 5       ]
  !  [      . . .   ]
  !  [        . . . ]
  subroutine Tpetra_CrsMatrix_CreateLabelledMat(comm,Mat)
    type(TeuchosComm), intent(in) :: comm
    type(TpetraCrsMatrix), intent(out) :: Mat
    type(TpetraMap) :: rowMap
    integer :: lclRow
    integer :: numProcs
    integer(global_ordinal_type) gblNumRows, gblNumCols, gblRow
    integer(global_ordinal_type), allocatable :: gblColInds(:)
    real(mag_type), allocatable :: vals(:)
    integer(global_ordinal_type), allocatable :: curGblColInds(:)
    real(mag_type), allocatable :: curVals(:)
    integer(size_type), parameter :: ithree=3

    ! create a Map
    numProcs = comm%getSize()
    gblNumRows = test_matrix_num_local() * numProcs
    rowMap = TpetraMap(gblNumRows, test_matrix_num_local(), comm)

    ! Leave room for three locally owned entries per row.
    ! We will only fill in two entries per row.
    Mat = TpetraCrsMatrix(rowMap, ithree)
    gblNumCols = gblNumRows

    allocate(gblColInds(2),vals(2))

    do lclRow = 1, test_matrix_num_local()
      gblRow = rowMap%getGlobalElement(lclRow)
      gblColInds(1) = 1 + modulo(gblRow + 0, gblNumCols)
      gblColInds(2) = 1 + modulo(gblRow + 1, gblNumCols)

      vals(1) = gblColInds(1)
      vals(2) = gblColInds(2)

      call Mat%insertGlobalValues(gblRow, gblColInds, vals)
    end do

    call rowMap%release()
    deallocate(gblColInds,vals)

    return
  end subroutine Tpetra_CrsMatrix_CreateLabelledMat

end module test_tpetra_crsmatrix_helper
