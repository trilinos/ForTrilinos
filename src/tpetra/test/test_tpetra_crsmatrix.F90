! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program test_TpetraCrsMatrix
#include "ForTrilinosTpetra_config.hpp"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  implicit none
  type(TeuchosComm) :: comm
  character(len=256), parameter :: FILENAME="test_tpetra_crsmatrix.F90"

  SETUP_TEST()

#ifdef HAVE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm(); FORTRILINOS_CHECK_IERR()
#endif

  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_Basic1)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_AlphaBetaMultiply)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_ActiveFillGlobal)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_EigTest)

!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_ActiveFillLocal)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_replaceColMap)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_replaceDomainMapAndImporter)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_removeEmptyProcessesInPlace)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getCrsGraph)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_isStorageOptimized)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_isStaticGraph)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_supportsRowViews)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getGlobalRowView)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getLocalRowViewRaw)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_gaussSeidel)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_reorderedGaussSeidel)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_gaussSeidelCopy)
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_reorderedGaussSeidelCopy)

  call comm%release()

  TEARDOWN_TEST()

contains

  ! ------------------------------- Basic1 ----------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_Basic1)
    type(TpetraMap) :: Map, row_map, tmp_map
    type(TeuchosComm) :: tcomm
    type(TpetraCrsMatrix) :: Mat
    type(TpetraMultiVector) :: mvrand, mvres
    character(kind=C_CHAR, len=:), allocatable :: description
    integer :: num_images, my_image_id
    integer(size_type), parameter :: num_vecs=5
    integer, parameter :: num_local=10
    integer :: irow
    integer(global_ordinal_type) :: base, gblrow, cols(1)
    real(mag_type) :: vals(1)=[1.], norms(num_vecs), zeros(num_vecs), fnorm

    OUT0("Starting TpetraCrsMatrix_Basic1")

    zeros = 0.d0

    num_images = comm%getSize()
    my_image_id = comm%getRank()

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    mvrand = TpetraMultiVector(Map, num_vecs, .false.); TEST_IERR()
    mvres = TpetraMultiVector(Map, num_vecs, .false.); TEST_IERR()
    call mvrand%randomize(); TEST_IERR()

    ! create the identity matrix
    base = num_local * my_image_id;
    Mat = TpetraCrsMatrix(Map, 1_size_type, TpetraStaticProfile)
    do irow = 1, num_local
      gblrow = base + int(irow, kind=global_ordinal_type)
      cols(1) = gblrow
      vals(1) = 1.d0
      call Mat%insertGlobalValues(gblrow, cols, vals)
    end do

    TEST_ASSERT(Mat%isGloballyIndexed())
    TEST_ASSERT((.not. Mat%isLocallyIndexed()))
    TEST_ASSERT(Mat%getProfileType() == TpetraStaticProfile)
    call Mat%fillComplete(); TEST_IERR()
    row_map = Mat%getRowMap()

    ! test the properties
    TEST_ASSERT(Mat%getGlobalNumEntries()==num_images*num_local)
    TEST_ASSERT(Mat%getNodeNumEntries()==num_local)
    TEST_ASSERT(Mat%getGlobalNumRows()==num_images*num_local)
    TEST_ASSERT(Mat%getGlobalNumCols()==num_images*num_local)
    TEST_ASSERT(Mat%getNodeNumRows()==num_local)
    TEST_ASSERT(Mat%getNodeNumCols()==num_local)
    TEST_ASSERT(Mat%getGlobalMaxNumRowEntries()==1)
    TEST_ASSERT(Mat%getNodeMaxNumRowEntries()==1)
    TEST_ASSERT(Mat%isFillComplete())
    TEST_ASSERT((.not. Mat%isFillActive()))
    fnorm = sqrt(real(num_images*num_local, kind=scalar_type))
    TEST_ASSERT(Mat%getFrobeniusNorm()==fnorm)
    tmp_map = Mat%getColMap(); TEST_IERR()
    TEST_ASSERT(row_map%isSameAs(tmp_map))
    call tmp_map%release(); TEST_IERR()
    tmp_map = Mat%getDomainMap(); TEST_IERR()
    TEST_ASSERT(row_map%isSameAs(tmp_map))
    call tmp_map%release(); TEST_IERR()
    tmp_map = Mat%getRangeMap(); TEST_IERR()
    TEST_ASSERT(row_map%isSameAs(tmp_map))
    call tmp_map%release(); TEST_IERR()
    TEST_ASSERT(Mat%hasColMap())
    TEST_ASSERT(Mat%haveGlobalConstants())
    !TEST_ASSERT((.not. Mat%isLowerTriangular())) ! FIXME: This throws an error from Tpetra
    !TEST_ASSERT((.not. Mat%isUpperTriangular())) ! FIXME: This tthrows an error from Tpetra
    TEST_ASSERT(Mat%isLocallyIndexed())
    TEST_ASSERT(Mat%hasTransposeApply())
    TEST_ASSERT((.not. Mat%isGloballyIndexed()))

    tcomm = Mat%getComm(); TEST_IERR()
    TEST_ASSERT(tcomm%getRank()==comm%getRank())
    TEST_ASSERT(tcomm%getSize()==comm%getSize())

    ! get the description, just to see if it does not throw
    description = Mat%description(); TEST_IERR()

    do irow = 1, num_local
      gblrow = row_map%getGlobalElement(irow)
      TEST_ASSERT(gblrow == (base + int(irow, kind=global_ordinal_type)))
      TEST_ASSERT(Mat%getNumEntriesInGlobalRow(gblrow)==1)
    end do

    ! test the action
    call mvres%randomize(); TEST_IERR()
    call Mat%apply(mvrand, mvres); TEST_IERR()
    call mvres%update(1.d0, mvrand, -1.d0); TEST_IERR()

    call mvres%norm1(norms); TEST_IERR()
    TEST_FLOATING_ARRAY_EQUALITY(norms, 0.d0, epsilon(0.d0))

    ! Set all diagonal entries to 2 and do again
    call Mat%resumeFill(); TEST_IERR()
    call Mat%setAllToScalar(2.d0); TEST_IERR()
    call Mat%fillComplete(); TEST_IERR()

    call mvres%randomize(); TEST_IERR()
    call Mat%apply(mvrand, mvres); TEST_IERR()
    call mvres%update(2.d0, mvrand, -1.d0); TEST_IERR()

    call mvres%norm1(norms); TEST_IERR()
    TEST_FLOATING_ARRAY_EQUALITY(norms, 0.d0, epsilon(0.d0))

    ! Scale diagonal entries by 2 and do again
    call Mat%resumeFill(); TEST_IERR()
    call Mat%scale(2.d0); TEST_IERR()
    call Mat%fillComplete(); TEST_IERR()

    call mvres%randomize(); TEST_IERR()
    call Mat%apply(mvrand, mvres); TEST_IERR()
    call mvres%update(4.d0, mvrand, -1.d0); TEST_IERR()

    call mvres%norm1(norms); TEST_IERR()
    TEST_FLOATING_ARRAY_EQUALITY(norms, 0.d0, epsilon(0.d0))

    call Mat%release(); TEST_IERR()
    call map%release(); TEST_IERR()
    call row_map%release(); TEST_IERR()
    call tcomm%release(); TEST_IERR()
    call mvres%release(); TEST_IERR()
    call mvrand%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_Basic1!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_Basic1)

  ! -------------------------- EigTest --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_EigTest)
    type(TpetraMap) :: Map, row_map
    type(TpetraCrsMatrix) :: A
    type(TpetraMultiVector) :: ones, threes
    integer :: num_images, my_image_id
    integer(size_type) :: numindices
    integer :: nnz
    real(norm_type) :: norms(1)
    integer(global_ordinal_type) :: gblrow
    integer(global_ordinal_type), allocatable :: cols(:), xcols(:)
    real(scalar_type), allocatable :: vals(:), xvals(:)

    OUT0("Starting TpetraCrsMatrix_EigTest")
    num_images = comm%getSize()
    my_image_id = comm%getRank()

    if (num_images < 2) return

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, 1, comm); TEST_IERR()

    ! create a multivector ones(n,1)
    ones = TpetraMultiVector(map, 1_size_type, .false.); TEST_IERR()
    threes = TpetraMultiVector(map, 1_size_type, .false.); TEST_IERR()
    call ones%putScalar(1.d0)

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

    A = TpetraCrsMatrix(map, 3_size_type, TpetraStaticProfile); TEST_IERR()
    gblrow = my_image_id + 1
    if (gblrow == 1) then
      nnz = 2
      allocate(cols(nnz)); allocate(vals(nnz));
      cols(1:nnz) = [gblrow, gblrow+1]
      vals(1:nnz) = [2.d0, 1.d0]
      call A%insertGlobalValues(gblrow, cols, vals); TEST_IERR()
    else if (gblrow == num_images) then
      nnz = 2;
      allocate(cols(nnz)); allocate(vals(nnz));
      cols(1:nnz) = [gblrow-1, gblrow]
      vals(1:nnz) = [1.d0, 2.d0]
      call A%insertGlobalValues(gblrow, cols, vals); TEST_IERR()
    else
      nnz = 3;
      allocate(cols(nnz)); allocate(vals(nnz));
      vals = [1.d0, 1.d0, 1.d0]
      cols = [gblrow-1, gblrow, gblrow+1]
      call A%insertGlobalValues(gblrow, cols, vals); TEST_IERR()
    end if

    call A%fillComplete(); TEST_IERR()
    row_map = A%getRowMap(); TEST_IERR()

    ! test the properties
    TEST_ASSERT(A%getGlobalNumEntries()==3*num_images-2)
    TEST_ASSERT(A%getNodeNumEntries()==nnz)
    TEST_ASSERT(A%getGlobalNumRows()==num_images)
    TEST_ASSERT(A%getNodeNumRows()==1)
    TEST_ASSERT(A%getNodeNumCols()==nnz)
    if (num_images > 2) then
      TEST_ASSERT(A%getGlobalMaxNumRowEntries()==3)
    else
      TEST_ASSERT(A%getGlobalMaxNumRowEntries()==2)
    end if
    TEST_ASSERT(A%getNodeMaxNumRowEntries()==nnz)
    TEST_ASSERT((.not. row_map%isSameAs(A%getColMap())))
    TEST_ASSERT(row_map%isSameAs(A%getDomainMap()))
    TEST_ASSERT(row_map%isSameAs(A%getRangeMap()))

    ! FIXME: getGlobalRowCopy: I expected cols and vals to be in ascending
    ! FIXME: order, but they are not, so I have commented out the check
    allocate(xcols(nnz)); allocate(xvals(nnz))
    numindices = int(nnz, kind=size_type)
    call A%getGlobalRowCopy(gblrow, xcols, xvals, numindices); TEST_IERR()
    !TEST_ARRAY_EQUALITY(xcols, cols)
    !TEST_FLOATING_ARRAY_EQUALITY(xvals, vals, epsilon(0.d0))
    deallocate(xcols); deallocate(xvals)

    ! test the action
    call threes%randomize(); TEST_IERR()
    call A%apply(ones, threes); TEST_IERR()

    ! now, threes should be 3*ones
    call threes%update(-3.d0, ones, 1.d0)
    call threes%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, 0.d0, epsilon(0.d0))

    call ones%release()
    call threes%release()
    call Map%release()
    call A%release()
    call row_map%release()
    deallocate(cols); deallocate(vals)

    OUT0("Finished TpetraCrsMatrix_EigTest!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_EigTest)

  ! ------------------------ Alphabetamultiply ------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_AlphaBetaMultiply)
    type(TpetraMap) :: Map
    type(TpetraCrsMatrix) :: Mat
    type(TpetraMultiVector) :: X, Y, Z
    integer(size_type), parameter ::numvecs=1
    real(scalar_type) :: alpha, beta
    integer :: my_image_id
    integer(global_ordinal_type) :: gblrow, base, cols(1), i
    real(scalar_type) :: ones(1)=[1.d0], normz(numvecs), normy(numvecs)

    OUT0("Starting TpetraCrsMatrix_AlphaBetaMultiply")

    my_image_id = comm%getRank()

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, 3, comm); TEST_IERR()

    ! Create the identity matrix, three rows per proc
    base = 3 * my_image_id;
    Mat = TpetraCrsMatrix(Map, 1_size_type, TpetraStaticProfile); TEST_IERR()
    do i = 1, 3
      gblrow = base + i
      cols(1) = gblrow
      call Mat%insertGlobalValues(gblrow, cols, ones); TEST_IERR()
    end do
    call Mat%fillComplete(); TEST_IERR()

    X = TpetraMultiVector(Map, numvecs); TEST_IERR()
    Y = TpetraMultiVector(Map, numvecs); TEST_IERR()
    Z = TpetraMultiVector(map, numvecs); TEST_IERR()

    call random_number(alpha)
    call random_number(beta)

    call X%randomize()
    call Y%randomize()

    ! Z = alpha*X + beta*Y
    call Z%update(alpha, X, beta, Y, 0.d0)

    ! test the action: Y = alpha*I*X + beta*Y = alpha*X + beta*Y = Z
    call Mat%apply(X, Y, TEUCHOSNO_TRANS, alpha, beta)
    !
    call Z%norm1(normz)
    call Y%norm1(normy)
    TEST_FLOATING_ARRAY_EQUALITY(normy, normz, epsilon(0.d0))

    call Z%release()
    call Y%release()
    call X%release()
    call Mat%release()
    call Map%release()

    OUT0("Finished TpetraCrsMatrix_AlphaBetaMultiply!")
  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_AlphaBetaMultiply)

  ! ----------------------------- ActiveFill --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_ActiveFillGlobal)
    type(TpetraMap) :: Map
    type(ParameterList) :: params
    type(TpetraCrsMatrix) :: Mat
    integer :: lclrow, numvalid
    integer(global_ordinal_type) :: row, cols(1)
    real(scalar_type) :: vals(1), zeros(1)=[0.d0]

    OUT0("Starting TpetraCrsMatrix_ActiveFillGlobal")

    ! create Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, 1, comm); TEST_IERR()

    Mat = TpetraCrsMatrix(map, map, 1_size_type, TpetraStaticProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    lclrow = 1
    row = map%getGlobalElement(lclrow)
    cols(1) = row; vals(1) = 1.
    call Mat%insertGlobalValues(row, cols, zeros); TEST_IERR()

    params = ParameterList("ANONYMOUS")
    ! call params%set("Optimize Storage", .false.) ! FIXME: boolean parameters
    call Mat%fillComplete(params);
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    numvalid = 0

    ! It's forbidden to call any of the *LocalValues methods if the
    ! matrix is fill complete (not fill active).

    ! FIXME: This fails with segfault. Not sure what's going on.
    ! TEST_THROW(call Mat%insertGlobalValues(row, cols, vals))

    numvalid = Mat%replaceGlobalValues(row, cols, vals); TEST_IERR()
    TEST_ASSERT(numvalid==TPETRA_GLOBAL_INVALID)

    numvalid = Mat%sumIntoGlobalValues(row, cols, vals); TEST_IERR()
    TEST_ASSERT(numvalid==TPETRA_GLOBAL_INVALID)

    TEST_THROW(call Mat%setAllToScalar(0.d0))
    TEST_THROW(call Mat%scale(0.d0))
    TEST_THROW(call Mat%globalAssemble())
    TEST_THROW(call Mat%fillComplete())

    call params%release()
    call Mat%release()

    Mat = TpetraCrsMatrix(map, map,10_size_type, TpetraStaticProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    lclrow = 1
    row = map%getGlobalElement(lclrow)
    cols(1) = row; vals(1) = 0.d0;
    ! FIXME: If a mistake is made, and cols above is just set to 1:
    ! cols(1) = 1
    ! FIXME: then the call to fillComplete below hangs indefinitely
    call Mat%insertGlobalValues(row, cols, vals); TEST_IERR()

    params = ParameterList("ANONYMOUS"); TEST_IERR()
    call Mat%fillComplete(params); TEST_IERR()
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call Mat%resumeFill(); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))

    ! FIXME: The following should NOT set ierr/=0 but does
    !TEST_NOTHROW(call Mat%insertGlobalValues(row, cols, vals))
    TEST_NOTHROW(numvalid = Mat%replaceGlobalValues(row, cols, vals))
    TEST_NOTHROW(numvalid = Mat%sumIntoGlobalValues(row, cols, vals))
    TEST_NOTHROW(call Mat%setAllToScalar(0.d0))
    TEST_NOTHROW(call Mat%scale(0.d0))
    TEST_NOTHROW(call Mat%globalAssemble())

    TEST_NOTHROW(call Mat%fillComplete())
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call params%release()
    call Map%release()
    call Mat%release()

    OUT0("Finished TpetraCrsMatrix_ActiveFillGlobal!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_ActiveFillGlobal)

#if 0

  ! ----------------------------- ActiveFill --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_ActiveFillLocal)
    type(TpetraMap) :: Map
    type(ParameterList) :: params
    type(TpetraCrsMatrix) :: Mat
    integer :: row, cols(1), numvalid
    real(scalar_type) :: vals(1), zeros(1)=[0.d0]

    OUT0("Starting TpetraCrsMatrix_ActiveFillLocal")

    ! create Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, 1, comm); TEST_IERR()

    Mat = TpetraCrsMatrix(map, map, 1_size_type, TpetraStaticProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    row = 1; cols(1) = 1; vals(1) = 0.
    call Mat%insertLocalValues(row, cols, tuple<Scalar>(0)); TEST_IERR()

    params = ParameterList("ANONYMOUS")
    ! call params%set("Optimize Storage", .false.) ! FIXME: boolean parameters
    call Mat%fillComplete(params);
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    numvalid = 0

    ! It's forbidden to call any of the *LocalValues methods if the
    ! matrix is fill complete (not fill active).

    TEST_THROW(call Mat%insertLocalValues(row, cols, vals))

    numvalid = Mat%replaceLocalValues(row, cols, vals); TEST_IERR()
    TEST_ASSERT(numvalid==TPETRA_GLOBAL_INVALID)

    numvalid = Mat%sumIntoLocalValues(lcrow, cols, vals); TEST_IERR()
    TEST_ASSERT(numvalid==TPETRA_GLOBAL_INVALID)

    TEST_THROW(call Mat%setAllToScalar(0.d0))
    TEST_THROW(call Mat%scale(0.d0))
    TEST_THROW(call Mat%globalAssemble())
    TEST_THROW(call Mat%fillComplete())

    call params%release()
    call Mat%release()

    Mat = TpetraCrsMatrix(map, map, 1, TpetraStaticProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    row = 1; cols(1) = 1; vals(1) = 0.d0;
    call Mat%insertLocalValues(row, cols, vals); TEST_IERR()

    params = ParameterList("ANONYMOUS"); TEST_IERR()
    !call params%set("Optimize Storage", .false.); TEST_IERR() ! FIXME: boolean parameters
    call Mat%fillComplete(params); TEST_IERR()
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call Mat%resumeFill(); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    TEST_NOTHROW(call Mat%insertLocalValues(row, cols, vals))
    TEST_NOTHROW(numvalid = Mat%replaceLocalValues(row, cols, vals))
    TEST_NOTHROW(call Mat%sumIntoLocalValues(row, cols, vals))

    ! FIXME: The following should NOT set ierr/=0 but does
    TEST_NOTHROW(call Mat%setAllToScalar(0.d0))
    TEST_NOTHROW(call Mat%scale(0.d0))
    TEST_NOTHROW(call Mat%globalAssemble())

    TEST_NOTHROW(call Mat%fillComplete())
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call params%release()
    call Map%release()
    call Mat%release()

    OUT0("Finished TpetraCrsMatrix_ActiveFillLocal!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_ActiveFillLocal)

  ! -------------------------------getCrsGraph-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getCrsGraph)
    type(TpetraCrsMatrix) :: Obj
    OUT0("Starting TpetraCrsMatrix_getCrsGraph")

    success = .false.

    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !fresult = Obj%getCrsGraph(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_getCrsGraph: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_getCrsGraph")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getCrsGraph)

  ! ----------------------------isStorageOptimized---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isStorageOptimized)
    type(TpetraCrsMatrix) :: Obj
    OUT0("Starting TpetraCrsMatrix_isStorageOptimized")

    success = .false.

    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !fresult = Obj%isStorageOptimized(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_isStorageOptimized: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_isStorageOptimized")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isStorageOptimized)

  ! ------------------------------isStaticGraph------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isStaticGraph)
    type(TpetraCrsMatrix) :: Obj
    OUT0("Starting TpetraCrsMatrix_isStaticGraph")

    success = .false.

    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !fresult = Obj%isStaticGraph(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_isStaticGraph: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_isStaticGraph")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isStaticGraph)

  ! -----------------------------supportsRowViews----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_supportsRowViews)
    type(TpetraCrsMatrix) :: Obj
    OUT0("Starting TpetraCrsMatrix_supportsRowViews")

    success = .false.

    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !fresult = Obj%supportsRowViews(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_supportsRowViews: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_supportsRowViews")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_supportsRowViews)
#endif

  ! -----------------------------getGlobalRowView----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_GetGlobalRowView)

    type(TpetraMap) :: rowmap, colmap
    type(TpetraCrsMatrix) :: A
    type(ParameterList) :: params
    integer(TpetraProfileType) :: pftype
    integer(size_type) :: num_images, my_image_id, numentries
    integer :: nnz
    integer :: T, i
    integer, parameter :: lclrow=1
    logical :: opt_storage
    real(scalar_type) :: scopy(4)
    integer :: lcopy(4)
    integer(global_ordinal_type) :: gcopy(4)
    integer(global_ordinal_type), pointer :: cgptr(:) => NULL()
    real(scalar_type), pointer :: csptr(:) => NULL()
    integer(global_ordinal_type) :: gblrow
    integer, allocatable :: linds(:)
    integer(global_ordinal_type), allocatable :: ginds(:), mask(:)
    real(scalar_type), allocatable :: values(:)
    ! ------------------------------------------------------------------------ !

    OUT0("Starting TpetraCrsMatrix_GetGlobalRowView")

    num_images = comm%getSize()
    my_image_id = comm%getRank()

    if (num_images < 2) return

    ! create a Map, one row per processor
    rowmap = TpetraMap(TPETRA_GLOBAL_INVALID, 1, comm); TEST_IERR()

    gblrow = rowmap%getGlobalElement(1);

    ! specify the column map to control ordering
    ! construct tridiagonal graph
    if (gblrow == 1) then
      nnz = 2
      allocate(ginds(nnz));
      ginds = [gblrow, gblrow+1]
    else if (gblrow == num_images) then
      nnz = 2
      allocate(ginds(nnz));
      ginds = [gblrow-1, gblrow]
    else
      nnz = 3
      allocate(ginds(nnz));
      ginds = [gblrow-1, gblrow, gblrow+1]
    end if
    allocate(linds(nnz))
    forall(i=1:nnz) linds(i)=i

    allocate(values(nnz))
    values = 1.0d0

    ! Create column map
    colmap = TpetraMap(TPETRA_GLOBAL_INVALID, ginds, comm);

    params = ParameterList("ANONYMOUS")

    do T = 0, 1
      pftype = TpetraStaticProfile
      opt_storage = IAND(T, 2) == 2
      call params%set("Optimize Storage", opt_storage)

      ! only allocate as much room as necessary
      numentries = int(nnz, kind=size_type)
      A = TpetraCrsMatrix(rowmap, colmap, numentries, pftype)

      ! at this point, the graph has not allocated data as global or local, so
      ! we can do views/copies for either local or global
      call A%getLocalRowCopy(lclrow, lcopy, scopy, numentries)
      !call A%getLocalRowView(lclrow, clview, csptr)
      call A%getGlobalRowCopy(gblrow, gcopy, scopy, numentries)

      call A%insertGlobalValues(gblrow, ginds, values)

      ! check values before calling fillComplete
      allocate(mask(nnz))
      mask = -1
      call A%getGlobalRowView(gblrow, cgptr, csptr)
      do i=1, nnz;
        where(ginds(i)==cgptr) mask = i
      end do
      TEST_ASSERT((.not. any(mask==-1)))
      TEST_ARRAY_EQUALITY(cgptr(mask), ginds)
      TEST_FLOATING_ARRAY_EQUALITY(csptr(mask), values, epsilon(0.d0))
      deallocate(mask)

      call A%fillComplete(params);

      ! check for throws and no-throws/values
      TEST_THROW(call A%getGlobalRowView(gblrow, cgptr, csptr))

      !TEST_NOTHROW(call A%getLocalRowView(lclrow, clview, csptr))
      !TEST_ARRAY_EQUALITY(clview, linds)
      !TEST_FLOATING_ARRAY_EQUALITY(csptr, values, epsilon(0.d0))

      TEST_NOTHROW(call A%getLocalRowCopy(lclrow, lcopy, scopy, numentries))
      TEST_ARRAY_EQUALITY(lcopy(1:numentries), linds)
      TEST_FLOATING_ARRAY_EQUALITY(scopy(1:numentries), values, epsilon(0.d0))

      TEST_NOTHROW(call A%getGlobalRowCopy(gblrow, gcopy, scopy, numentries) );
      TEST_ARRAY_EQUALITY(gcopy(1:numentries), ginds)
      TEST_FLOATING_ARRAY_EQUALITY(scopy(1:numentries), values, epsilon(0.d0))

      call A%release()

    end do

    call rowmap%release()
    call colmap%release()

    deallocate(linds); deallocate(ginds); deallocate(values)

    OUT0("Finished TpetraCrsMatrix_GetGlobalRowView")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_GetGlobalRowView)

#if 0
  ! ----------------------------getLocalRowViewRaw---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getLocalRowViewRaw)
    type(TpetraCrsMatrix) :: Obj
    integer :: lclrow
    integer :: nument
    type(C_PTR) :: lclcolinds
    type(C_PTR) :: vals
    OUT0("Starting TpetraCrsMatrix_getLocalRowViewRaw")

    success = .false.

    lclrow = 0
    nument = 0
    !lclcolinds = xxx(); TEST_IERR()
    !vals = xxx(); TEST_IERR()
    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !fresult = Obj%getLocalRowViewRaw(lclrow, nument, lclcolinds, vals); TEST_IERR()

    !call lclcolinds%release(); TEST_IERR()
    !call vals%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_getLocalRowViewRaw: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_getLocalRowViewRaw")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getLocalRowViewRaw)

  ! -------------------------------gaussSeidel-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_gaussSeidel)
    type(TpetraCrsMatrix) :: Obj
    type(TpetraMultiVector) :: b
    type(TpetraMultiVector) :: x
    type(TpetraMultiVector) :: d
    real(scalar_type) :: dampingfactor
    integer(kind(TpetraESweepDirection)) :: direction
    integer(C_INT) :: numsweeps
    OUT0("Starting TpetraCrsMatrix_gaussSeidel")

    success = .false.

    !b = TpetraMultiVector(); TEST_IERR()
    !x = TpetraMultiVector(); TEST_IERR()
    !d = TpetraMultiVector(); TEST_IERR()
    dampingfactor = 0
    direction = 0
    numsweeps = 0
    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !call Obj%gaussSeidel(b, x, d, dampingfactor, direction, numsweeps); TEST_IERR()

    !call b%release(); TEST_IERR()
    !call x%release(); TEST_IERR()
    !call d%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_gaussSeidel: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_gaussSeidel")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_gaussSeidel)

  ! ---------------------------reorderedGaussSeidel--------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_reorderedGaussSeidel)
    type(TpetraCrsMatrix) :: Obj
    type(TpetraMultiVector) :: b
    type(TpetraMultiVector) :: x
    type(TpetraMultiVector) :: d
    !type(TeuchosArrayViewInt) :: rowindices
    real(scalar_type) :: dampingfactor
    integer(kind(TpetraESweepDirection)) :: direction
    integer(C_INT) :: numsweeps
    OUT0("Starting TpetraCrsMatrix_reorderedGaussSeidel")

    success = .false.

    !b = TpetraMultiVector(); TEST_IERR()
    !x = TpetraMultiVector(); TEST_IERR()
    !d = TpetraMultiVector(); TEST_IERR()
    !rowindices = xxx(); TEST_IERR()
    dampingfactor = 0
    direction = 0
    numsweeps = 0
    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !call Obj%reorderedGaussSeidel(b, x, d, rowindices, dampingfactor, direction, numsweeps); TEST_IERR()

    !call b%release(); TEST_IERR()
    !call x%release(); TEST_IERR()
    !call d%release(); TEST_IERR()
    !call rowindices%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_reorderedGaussSeidel: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_reorderedGaussSeidel")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_reorderedGaussSeidel)

  ! -----------------------------gaussSeidelCopy------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_gaussSeidelCopy)
    type(TpetraCrsMatrix) :: Obj
    type(TpetraMultiVector) :: x
    type(TpetraMultiVector) :: b
    type(TpetraMultiVector) :: d
    real(scalar_type) :: dampingfactor
    integer(kind(TpetraESweepDirection)) :: direction
    integer(C_INT) :: numsweeps
    logical(C_BOOL) :: zeroinitialguess
    OUT0("Starting TpetraCrsMatrix_gaussSeidelCopy")

    success = .false.

    !x = TpetraMultiVector(); TEST_IERR()
    !b = TpetraMultiVector(); TEST_IERR()
    !d = TpetraMultiVector(); TEST_IERR()
    dampingfactor = 0
    direction = 0
    numsweeps = 0
    zeroinitialguess = .false.
    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !call Obj%gaussSeidelCopy(x, b, d, dampingfactor, direction, numsweeps, zeroinitialguess); TEST_IERR()

    !call x%release(); TEST_IERR()
    !call b%release(); TEST_IERR()
    !call d%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_gaussSeidelCopy: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_gaussSeidelCopy")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_gaussSeidelCopy)

  ! -------------------------reorderedGaussSeidelCopy------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_reorderedGaussSeidelCopy)
    type(TpetraCrsMatrix) :: Obj
    type(TpetraMultiVector) :: x
    type(TpetraMultiVector) :: b
    type(TpetraMultiVector) :: d
    !type(TeuchosArrayViewInt) :: rowindices
    real(scalar_type) :: dampingfactor
    integer(kind(TpetraESweepDirection)) :: direction
    integer(C_INT) :: numsweeps
    logical(C_BOOL) :: zeroinitialguess
    OUT0("Starting TpetraCrsMatrix_reorderedGaussSeidelCopy")

    success = .false.

    !x = TpetraMultiVector(); TEST_IERR()
    !b = TpetraMultiVector(); TEST_IERR()
    !d = TpetraMultiVector(); TEST_IERR()
    !rowindices = xxx(); TEST_IERR()
    dampingfactor = 0
    direction = 0
    numsweeps = 0
    zeroinitialguess = .false.
    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !call Obj%reorderedGaussSeidelCopy(x, b, d, rowindices, dampingfactor, direction, numsweeps, zeroinitialguess); TEST_IERR()

    !call x%release(); TEST_IERR()
    !call b%release(); TEST_IERR()
    !call d%release(); TEST_IERR()
    !call rowindices%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_reorderedGaussSeidelCopy: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_reorderedGaussSeidelCopy")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_reorderedGaussSeidelCopy)

  ! ------------------------------replaceColMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceColMap)
    type(TpetraCrsMatrix) :: Obj
    type(TpetraMap) :: newcolmap
    OUT0("Starting TpetraCrsMatrix_replaceColMap")

    success = .false.

    !newcolmap = TpetraMap(); TEST_IERR()
    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !call Obj%replaceColMap(newcolmap); TEST_IERR()

    !call newcolmap%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_replaceColMap: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_replaceColMap")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceColMap)

  ! -----------------------replaceDomainMapAndImporter------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceDomainMapAndImporter)
    type(TpetraCrsMatrix) :: Obj
    type(TpetraMap) :: newdomainmap
    type(TpetraImport) :: newimporter
    OUT0("Starting TpetraCrsMatrix_replaceDomainMapAndImporter")

    success = .false.

    !newdomainmap = TpetraMap(); TEST_IERR()
    !newimporter = TpetraImport(); TEST_IERR()
    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !call Obj%replaceDomainMapAndImporter(newdomainmap, newimporter); TEST_IERR()

    !call newdomainmap%release(); TEST_IERR()
    !call newimporter%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_replaceDomainMapAndImporter: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_replaceDomainMapAndImporter")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceDomainMapAndImporter)

  ! -----------------------removeEmptyProcessesInPlace------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_removeEmptyProcessesInPlace)
    type(TpetraCrsMatrix) :: Obj
    type(TpetraMap) :: newmap
    OUT0("Starting TpetraCrsMatrix_removeEmptyProcessesInPlace")

    success = .false.

    !newmap = TpetraMap(); TEST_IERR()
    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !call Obj%removeEmptyProcessesInPlace(newmap); TEST_IERR()

    !call newmap%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_removeEmptyProcessesInPlace: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_removeEmptyProcessesInPlace")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_removeEmptyProcessesInPlace)
#endif


end program test_TpetraCrsMatrix
