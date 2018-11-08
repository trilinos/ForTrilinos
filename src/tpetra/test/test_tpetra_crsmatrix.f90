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
  character(len=256), parameter :: FILENAME="test_tpetra_crsmatrix.f90"

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
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getGlobalRowView)
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
    logical(c_bool), parameter :: false=.false., true=.true.
    integer(size_type) :: num_images, my_image_id
    integer(size_type), parameter :: num_local=10, num_vecs=5
    integer(size_type), parameter :: ione=1
    integer(local_ordinal_type) :: irow
    integer(global_ordinal_type) :: base, gblrow, cols(1)
    real(mag_type) :: vals(1)=[1.], norms(num_vecs), zeros(num_vecs), fnorm
    real(scalar_type), parameter :: zero=0., one=1., negone=-1., two=2., four=4.
    real(scalar_type), allocatable :: a1(:), a2(:)
    integer(size_type) :: lda

    OUT0("Starting TpetraCrsMatrix_Basic1")

    zeros = zero

    num_images = comm%getSize()
    my_image_id = comm%getRank()

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    mvrand = TpetraMultiVector(Map, num_vecs, false); TEST_IERR()
    mvres = TpetraMultiVector(Map, num_vecs, false); TEST_IERR()
    call mvrand%randomize(); TEST_IERR()

    ! create the identity matrix
    base = num_local * my_image_id;
    Mat = TpetraCrsMatrix(Map, ione, TpetraDynamicProfile)
    do irow = 1, num_local
      gblrow = base + int(irow, kind=global_ordinal_type)
      cols(1) = gblrow
      vals(1) = one
      call Mat%insertGlobalValues(gblrow, cols, vals)
    end do

    TEST_ASSERT(Mat%isGloballyIndexed())
    TEST_ASSERT((.not. Mat%isLocallyIndexed()))
    TEST_ASSERT(Mat%getProfileType() == TpetraDynamicProfile)
    call Mat%fillComplete(); TEST_IERR()
    row_map = Mat%getRowMap()

    ! test the properties
    TEST_ASSERT(Mat%getGlobalNumEntries()==num_images*num_local)
    TEST_ASSERT(Mat%getNodeNumEntries()==num_local)
    TEST_ASSERT(Mat%getGlobalNumRows()==num_images*num_local)
    TEST_ASSERT(Mat%getGlobalNumCols()==num_images*num_local)
    TEST_ASSERT(Mat%getNodeNumRows()==num_local)
    TEST_ASSERT(Mat%getNodeNumCols()==num_local)
    TEST_ASSERT(Mat%getGlobalNumDiags()==num_images*num_local)
    TEST_ASSERT(Mat%getNodeNumDiags()==num_local)
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
    call mvres%update(one, mvrand, negone); TEST_IERR()

    call mvres%norm1(norms); TEST_IERR()
    TEST_FLOATING_ARRAY_EQUALITY(norms, zero, epsilon(zero))

    ! Set all diagonal entries to 2 and do again
    call Mat%resumeFill(); TEST_IERR()
    call Mat%setAllToScalar(two); TEST_IERR()
    call Mat%fillComplete(); TEST_IERR()

    call mvres%randomize(); TEST_IERR()
    call Mat%apply(mvrand, mvres); TEST_IERR()
    call mvres%update(two, mvrand, negone); TEST_IERR()

    call mvres%norm1(norms); TEST_IERR()
    TEST_FLOATING_ARRAY_EQUALITY(norms, zero, epsilon(zero))

    ! Scale diagonal entries by 2 and do again
    call Mat%resumeFill(); TEST_IERR()
    call Mat%scale(two); TEST_IERR()
    call Mat%fillComplete(); TEST_IERR()

    call mvres%randomize(); TEST_IERR()
    call Mat%apply(mvrand, mvres); TEST_IERR()
    call mvres%update(four, mvrand, negone); TEST_IERR()

    call mvres%norm1(norms); TEST_IERR()
    TEST_FLOATING_ARRAY_EQUALITY(norms, zero, epsilon(zero))

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
    integer(size_type), parameter :: izero=0, ione=1, ithree=3
    integer(size_type) :: num_images, my_image_id, numindices
    integer(local_ordinal_type) :: nnz
    logical(c_bool), parameter :: false=.false., true=.true.
    real(scalar_type), parameter :: zero=0., one=1., two=2., negthree=-3.
    real(norm_type) :: norms(1)
    integer(global_ordinal_type) :: gblrow
    integer(global_ordinal_type), allocatable :: cols(:), xcols(:)
    real(scalar_type), allocatable :: vals(:), xvals(:)

    OUT0("Starting TpetraCrsMatrix_EigTest")
    num_images = comm%getSize()
    my_image_id = comm%getRank()

    if (num_images < 2) return

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, ione, comm); TEST_IERR()

    ! create a multivector ones(n,1)
    ones = TpetraMultiVector(map, ione, false); TEST_IERR()
    threes = TpetraMultiVector(map, ione, false); TEST_IERR()
    call ones%putScalar(one)

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

    A = TpetraCrsMatrix(map, izero, TpetraDynamicProfile); TEST_IERR()
    gblrow = my_image_id + 1
    if (gblrow == 1) then
      nnz = 2
      allocate(cols(nnz)); allocate(vals(nnz));
      cols(1:nnz) = [gblrow, gblrow+1]
      vals(1:nnz) = [two, one]
      call A%insertGlobalValues(gblrow, cols, vals); TEST_IERR()
    else if (gblrow == num_images) then
      nnz = 2;
      allocate(cols(nnz)); allocate(vals(nnz));
      cols(1:nnz) = [gblrow-1, gblrow]
      vals(1:nnz) = [one, two]
      call A%insertGlobalValues(gblrow, cols, vals); TEST_IERR()
    else
      nnz = 3;
      allocate(cols(nnz)); allocate(vals(nnz));
      vals = [one, one, one]
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
    TEST_ASSERT(A%getGlobalNumDiags()==num_images)
    TEST_ASSERT(A%getNodeNumDiags()==1)
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
    !TEST_FLOATING_ARRAY_EQUALITY(xvals, vals, epsilon(zero))
    deallocate(xcols); deallocate(xvals)

    ! test the action
    call threes%randomize(); TEST_IERR()
    call A%apply(ones, threes); TEST_IERR()

    ! now, threes should be 3*ones
    call threes%update(negthree, ones, one)
    call threes%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, zero, epsilon(zero))

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
    integer(size_type), parameter :: ione=1, ithree=3, numvecs=1
    real(scalar_type), parameter :: zero=0., one=1.
    real(scalar_type) :: alpha, beta
    integer(size_type) :: my_image_id
    integer(global_ordinal_type) :: gblrow, base, cols(1), i
    real(scalar_type) :: ones(1)=[one], normz(numvecs), normy(numvecs)

    OUT0("Starting TpetraCrsMatrix_AlphaBetaMultiply")

    my_image_id = comm%getRank()

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, ithree, comm); TEST_IERR()

    ! Create the identity matrix, three rows per proc
    base = 3 * my_image_id;
    Mat = TpetraCrsMatrix(Map, ione, TpetraDynamicProfile); TEST_IERR()
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
    call Z%update(alpha, X, beta, Y, zero)

    ! test the action: Y = alpha*I*X + beta*Y = alpha*X + beta*Y = Z
    call Mat%apply(X, Y, TEUCHOSNO_TRANS, alpha, beta)
    !
    call Z%norm1(normz)
    call Y%norm1(normy)
    TEST_FLOATING_ARRAY_EQUALITY(normy, normz, epsilon(zero))

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
    integer(size_type), parameter :: izero=0, ione=1
    logical(c_bool), parameter :: true=.true., false=.false.
    integer(local_ordinal_type) :: lclrow, numvalid
    integer(global_ordinal_type) :: row, cols(1)
    real(scalar_type), parameter :: zero=0.
    real(scalar_type) :: vals(1), zeros(1)=[zero]

    OUT0("Starting TpetraCrsMatrix_ActiveFillGlobal")

    ! create Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, ione, comm); TEST_IERR()

    Mat = TpetraCrsMatrix(map, map, izero, TpetraDynamicProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    lclrow = 1
    row = map%getGlobalElement(lclrow)
    cols(1) = row; vals(1) = 0.
    call Mat%insertGlobalValues(row, cols, zeros); TEST_IERR()

    params = ParameterList("ANONOMOUS")
    ! call params%set("Optimize Storage", false) ! FIXME: boolean parameters
    call Mat%fillComplete(params);
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    numvalid = 0

    ! It's forbidden to call any of the *LocalValues methods if the
    ! matrix is fill complete (not fill active).

    TEST_THROW(call Mat%insertGlobalValues(row, cols, vals))

    numvalid = Mat%replaceGlobalValues(row, cols, vals); TEST_IERR()
    TEST_ASSERT(numvalid==TPETRA_GLOBAL_INVALID)

    numvalid = Mat%sumIntoGlobalValues(row, cols, vals); TEST_IERR()
    TEST_ASSERT(numvalid==TPETRA_GLOBAL_INVALID)

    TEST_THROW(call Mat%setAllToScalar(zero))
    TEST_THROW(call Mat%scale(zero))
    TEST_THROW(call Mat%globalAssemble())
    TEST_THROW(call Mat%fillComplete())

    call params%release()
    call Mat%release()

    Mat = TpetraCrsMatrix(map, map, izero, TpetraDynamicProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    lclrow = 1
    row = map%getGlobalElement(lclrow)
    cols(1) = row; vals(1) = zero;
    ! FIXME: If a mistake is made, and cols above is just set to 1:
    ! cols(1) = 1
    ! FIXME: then the call to fillComplete below hangs indefinitely
    call Mat%insertGlobalValues(row, cols, vals); TEST_IERR()

    params = ParameterList("ANONOMOUS"); TEST_IERR()
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
    TEST_NOTHROW(call Mat%setAllToScalar(zero))
    TEST_NOTHROW(call Mat%scale(zero))
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
    integer(size_type), parameter :: izero=0, ione=1
    logical(c_bool), parameter :: true=.true., false=.false.
    integer(local_ordinal_type) :: row, cols(1), numvalid
    real(scalar_type), parameter :: zero=0.
    real(scalar_type) :: vals(1), zeros(1)=[zero]

    OUT0("Starting TpetraCrsMatrix_ActiveFillLocal")

    ! create Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, ione, comm); TEST_IERR()

    Mat = TpetraCrsMatrix(map, map, izero, TpetraDynamicProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    row = 1; cols(1) = 1; vals(1) = 0.
    call Mat%insertLocalValues(row, cols, tuple<Scalar>(0)); TEST_IERR()

    params = ParameterList("ANONOMOUS")
    ! call params%set("Optimize Storage", false) ! FIXME: boolean parameters
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

    TEST_THROW(call Mat%setAllToScalar(zero))
    TEST_THROW(call Mat%scale(zero))
    TEST_THROW(call Mat%globalAssemble())
    TEST_THROW(call Mat%fillComplete())

    call params%release()
    call Mat%release()

    Mat = TpetraCrsMatrix(map, map, izero, TpetraDynamicProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    row = 1; cols(1) = 1; vals(1) = zero;
    call Mat%insertLocalValues(row, cols, vals); TEST_IERR()

    params = ParameterList("ANONOMOUS"); TEST_IERR()
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
    TEST_NOTHROW(call Mat%setAllToScalar(zero))
    TEST_NOTHROW(call Mat%scale(zero))
    TEST_NOTHROW(call Mat%globalAssemble())

    TEST_NOTHROW(call Mat%fillComplete())
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call params%release()
    call Map%release()
    call Mat%release()

    OUT0("Finished TpetraCrsMatrix_ActiveFillLocal!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_ActiveFillLocal)
#endif

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

  ! -----------------------------getGlobalRowView----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalRowView)
    type(TpetraCrsMatrix) :: Obj
    integer(global_ordinal_type) :: globalrow
    !type(TeuchosArrayViewLongLongConst) :: indices
    !type(TeuchosArrayViewDoubleConst) :: values
    OUT0("Starting TpetraCrsMatrix_getGlobalRowView")

    success = .false.

    globalrow = 0
    !indices = xxx(); TEST_IERR()
    !values = xxx(); TEST_IERR()
    !Obj = TpetraCrsMatrix(); TEST_IERR()
    !call Obj%getGlobalRowView(globalrow, indices, values); TEST_IERR()

    !call indices%release(); TEST_IERR()
    !call values%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_getGlobalRowView: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_getGlobalRowView")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalRowView)

  ! ----------------------------getLocalRowViewRaw---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getLocalRowViewRaw)
    type(TpetraCrsMatrix) :: Obj
    integer(local_ordinal_type) :: lclrow
    integer(local_ordinal_type) :: nument
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


end program test_TpetraCrsMatrix
