! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program test_TpetraCrsMatrix
#include "ForTrilinos_config.h"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra
  use test_tpetra_crsmatrix_helper

  implicit none
  type(TeuchosComm) :: comm
  character(len=256), parameter :: FILENAME="test_tpetra_crsmatrix.F90"
  integer, parameter :: dp = kind(0.d0)

  SETUP_TEST()

#if FORTRILINOS_USE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm(); FORTRILINOS_CHECK_IERR()
#endif

  !Fat tests
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_Basic1)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_AlphaBetaMultiply)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_ActiveFillGlobal)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_EigTest)

  !Unit tests consistent with scripts/autogen_tests.py
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getComm)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getRowMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getColMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getDomainMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getRangeMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getCrsGraph)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getGlobalNumRows)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getGlobalNumCols)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getNodeNumRows)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getNodeNumCols)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getGlobalNumEntries)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getNodeNumEntries)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getNumEntriesInGlobalRow)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getNumEntriesInLocalRow)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getGlobalMaxNumRowEntries)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getNodeMaxNumRowEntries)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getGlobalRowCopy)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getLocalRowCopy)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getGlobalRowView)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getProfileType)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getFrobeniusNorm)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getAllValues)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_setAllValues)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_insertGlobalValues)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_insertLocalValues)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_replaceGlobalValues)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_setAllToScalar)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_scale)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_globalAssemble)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_replaceColMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_replaceDomainMapAndImporter)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_removeEmptyProcessesInPlace)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_hasColMap)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_isLocallyIndexed)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_isGloballyIndexed)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_isFillComplete)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_isFillActive)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_isStorageOptimized)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_isStaticGraph)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_supportsRowViews)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_hasTransposeApply)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_description)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_computeGlobalConstants)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_haveGlobalConstants)

  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_gaussSeidel)
  !TODO-implement ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_gaussSeidelCopy)

  call comm%release();  TEST_IERR()

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
    integer :: irow
    integer(global_ordinal_type) :: base, gblrow, cols(1)
    real(mag_type) :: vals(1)=[1.], norms(num_vecs), zeros(num_vecs), fnorm

    OUT0("Starting TpetraCrsMatrix_Basic1")

    zeros = 0.d0

    num_images = comm%getSize()
    my_image_id = comm%getRank()

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, test_matrix_num_local(), comm); TEST_IERR()
    mvrand = TpetraMultiVector(Map, num_vecs, .false.); TEST_IERR()
    mvres = TpetraMultiVector(Map, num_vecs, .false.); TEST_IERR()
    call mvrand%randomize(); TEST_IERR()

    ! create the identity matrix
    base = test_matrix_num_local() * my_image_id;
    Mat = TpetraCrsMatrix(Map, 1_size_type, TpetraStaticProfile)
    do irow = 1, test_matrix_num_local()
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
    TEST_ASSERT(Mat%getGlobalNumEntries()==num_images*test_matrix_num_local())
    TEST_ASSERT(Mat%getNodeNumEntries()==test_matrix_num_local())
    TEST_ASSERT(Mat%getGlobalNumRows()==num_images*test_matrix_num_local())
    TEST_ASSERT(Mat%getGlobalNumCols()==num_images*test_matrix_num_local())
    TEST_ASSERT(Mat%getNodeNumRows()==test_matrix_num_local())
    TEST_ASSERT(Mat%getNodeNumCols()==test_matrix_num_local())
    TEST_ASSERT(Mat%getGlobalMaxNumRowEntries()==1)
    TEST_ASSERT(Mat%getNodeMaxNumRowEntries()==1)
    TEST_ASSERT(Mat%isFillComplete())
    TEST_ASSERT((.not. Mat%isFillActive()))
    fnorm = sqrt(real(num_images*test_matrix_num_local(), kind=scalar_type))
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
    TEST_ASSERT(Mat%isLocallyIndexed())
    TEST_ASSERT(Mat%hasTransposeApply())
    TEST_ASSERT((.not. Mat%isGloballyIndexed()))

    tcomm = Mat%getComm(); TEST_IERR()
    TEST_ASSERT(tcomm%getRank()==comm%getRank())
    TEST_ASSERT(tcomm%getSize()==comm%getSize())

    ! get the description, just to see if it does not throw
    description = Mat%description(); TEST_IERR()

    do irow = 1, test_matrix_num_local()
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
    type(TpetraMap) :: row_map, tmp_map
    type(TpetraCrsMatrix) :: A
    type(TpetraMultiVector) :: ones, threes
    integer :: num_images, my_image_id, nnz, nrow, nrm
    integer(size_type) :: numindices
    real(norm_type) :: norms(1)
    integer(global_ordinal_type) :: gblrow
    integer(global_ordinal_type), allocatable :: xcols(:)
    real(scalar_type), allocatable :: vals(:), xvals(:)

    OUT0("Starting TpetraCrsMatrix_EigTest")
    num_images = comm%getSize()
    my_image_id = comm%getRank()

    call TPetra_CrsMatrix_CreateTestMatrix_A(comm,A)
    call A%fillComplete(); TEST_IERR()
    row_map = A%getRowMap(); TEST_IERR()

    ! Must match TPetra_CrsMatrix_CreateTestMatrix
    nnz=3
    nrow=4
    if (my_image_id==0 .or. my_image_id==num_images-1) then
       nrm=1
    else
       nrm=0
    endif
    if (num_images==1) nrm=2

    ! test the properties
    TEST_ASSERT(A%getGlobalNumEntries()==nnz*nrow*num_images-2)
    TEST_ASSERT(A%getNodeNumEntries()==nnz*nrow-nrm)
    TEST_ASSERT(A%getGlobalNumRows()==num_images*nrow)
    TEST_ASSERT(A%getNodeNumRows()==nrow)
    TEST_ASSERT(A%getNodeNumCols()==nnz+(nrow-1)-nrm)
    TEST_ASSERT(A%getGlobalMaxNumRowEntries()==nnz)
    TEST_ASSERT(A%getNodeMaxNumRowEntries()==nnz)
    ! Use tmp_map so it can be released
    tmp_map=A%getColMap()
    if (num_images==1) then
      TEST_ASSERT(row_map%isSameAs(tmp_map))        ! Same in serial
    else
      TEST_ASSERT(.not. row_map%isSameAs(tmp_map))  ! not same on all procs
    endif
    tmp_map=A%getDomainMap()
    TEST_ASSERT(row_map%isSameAs(tmp_map))
    tmp_map=A%getRangeMap()
    TEST_ASSERT(row_map%isSameAs(tmp_map))

    ! create a multivector ones(n,1)
    ones = TpetraMultiVector(row_map, 1_size_type, .false.); TEST_IERR()
    threes = TpetraMultiVector(row_map, 1_size_type, .false.); TEST_IERR()
    call ones%putScalar(1.d0)
    ! test the action
    call threes%randomize(); TEST_IERR()
    call A%apply(ones, threes); TEST_IERR()

    ! now, threes should be 3*ones
    call threes%update(-3.d0, ones, 1.d0)
    call threes%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, 0.d0, epsilon(0.d0))

    call ones%release();  TEST_IERR()
    call threes%release();  TEST_IERR()
    call A%release();  TEST_IERR()
    call row_map%release();  TEST_IERR()
    call tmp_map%release();  TEST_IERR()

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

    call Z%release();  TEST_IERR()
    call Y%release();  TEST_IERR()
    call X%release();  TEST_IERR()
    call Mat%release();  TEST_IERR()
    call Map%release();  TEST_IERR()

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

    call params%release();  TEST_IERR()
    call Mat%release();  TEST_IERR()

    Mat = TpetraCrsMatrix(map, map,10_size_type, TpetraStaticProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    lclrow = 1
    row = map%getGlobalElement(lclrow)
    cols(1) = row; vals(1) = 0.d0;
    call Mat%insertGlobalValues(row, cols, vals); TEST_IERR()

    params = ParameterList("ANONYMOUS"); TEST_IERR()
    call Mat%fillComplete(params); TEST_IERR()
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call Mat%resumeFill(); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))

    TEST_NOTHROW(numvalid = Mat%replaceGlobalValues(row, cols, vals))
    TEST_NOTHROW(numvalid = Mat%sumIntoGlobalValues(row, cols, vals))
    TEST_NOTHROW(call Mat%setAllToScalar(0.d0))
    TEST_NOTHROW(call Mat%scale(0.d0))
    TEST_NOTHROW(call Mat%globalAssemble())

    TEST_NOTHROW(call Mat%fillComplete())
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call params%release();  TEST_IERR()
    call Map%release();  TEST_IERR()
    call Mat%release();  TEST_IERR()

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

    call params%release();  TEST_IERR()
    call Mat%release();  TEST_IERR()

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
    TEST_NOTHROW(call Mat%setAllToScalar(zero))
    TEST_NOTHROW(call Mat%scale(zero))
    TEST_NOTHROW(call Mat%globalAssemble())

    TEST_NOTHROW(call Mat%fillComplete())
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call params%release();  TEST_IERR()
    call Map%release();  TEST_IERR()
    call Mat%release();  TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_ActiveFillLocal!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_ActiveFillLocal)

#endif

  ! ---------------------------------getComm---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getComm)
    type(TpetraCrsMatrix) :: Mat
    type(TeuchosComm) :: tcomm

    OUT0("Starting TpetraCrsMatrix_getComm!")

    ! get identity matrix
    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    tcomm = Mat%getComm()

    TEST_ASSERT(tcomm%getRank() == comm%getRank())
    TEST_ASSERT(tcomm%getSize() == comm%getSize())

    call tcomm%release(); TEST_IERR()
    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getComm!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getComm)

  ! --------------------------------getRowMap--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getRowMap)
    type(TpetraMap) :: Map, tmpmap
    type(TpetraCrsMatrix) :: Mat
    integer :: lnum_local=2

    OUT0("Starting TpetraCrsMatrix_getRowMap!")

    ! create a Map
    Map = TpetraMap(TPETRA_GLOBAL_INVALID, lnum_local, comm); TEST_IERR()
    Mat = TpetraCrsMatrix(Map, 1_size_type, TpetraStaticProfile)

    tmpmap = Mat%getRowMap()
    TEST_EQUALITY(tmpmap%getNodeNumElements(), Map%getNodeNumElements())
    TEST_ASSERT(tmpmap%isSameAs(Map))

    call tmpmap%release(); TEST_IERR()
    call Map%release(); TEST_IERR()
    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getRowMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getRowMap)

  ! --------------------------------getColMap--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getColMap)
    type(TpetraMap) :: Map, tmpmap
    type(TpetraCrsMatrix) :: Mat
    integer :: lnum_local=2

    OUT0("Starting TpetraCrsMatrix_getColMap!")

    ! create a Map
    Map = TpetraMap(TPETRA_GLOBAL_INVALID, lnum_local, comm); TEST_IERR()
    Mat = TpetraCrsMatrix(Map, Map, 1_size_type, TpetraStaticProfile)
    TEST_NOTHROW(call Mat%fillComplete())

    tmpmap = Mat%getColMap()
    TEST_EQUALITY(tmpmap%getNodeNumElements(), Map%getNodeNumElements())

    call tmpmap%release(); TEST_IERR()
    call Map%release(); TEST_IERR()
    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getColMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getColMap)

  ! -------------------------------getDomainMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getDomainMap)
    type(TpetraMap) :: Map, tmpmap
    type(TpetraCrsMatrix) :: Mat
    integer :: lnum_local=2

    OUT0("Starting TpetraCrsMatrix_getDomainMap!")

    ! create a Map
    Map = TpetraMap(TPETRA_GLOBAL_INVALID, lnum_local, comm); TEST_IERR()
    Mat = TpetraCrsMatrix(Map, Map, 1_size_type, TpetraStaticProfile)
    TEST_NOTHROW(call Mat%fillComplete())

    tmpmap = Mat%getDomainMap(); TEST_IERR()
    TEST_ASSERT(Map%isSameAs(tmpmap))

    call tmpmap%release(); TEST_IERR()
    call Map%release(); TEST_IERR()
    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getDomainMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getDomainMap)

  ! -------------------------------getRangeMap-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getRangeMap)
    type(TpetraMap) :: Map, tmpmap
    type(TpetraCrsMatrix) :: Mat
    integer :: lnum_local=2

    OUT0("Starting TpetraCrsMatrix_getRangeMap!")

    ! create a Map
    Map = TpetraMap(TPETRA_GLOBAL_INVALID, lnum_local, comm); TEST_IERR()
    Mat = TpetraCrsMatrix(Map, Map, 1_size_type, TpetraStaticProfile)
    TEST_NOTHROW(call Mat%fillComplete())

    tmpmap = Mat%getRangeMap(); TEST_IERR()
    TEST_ASSERT(Map%isSameAs(tmpmap))

    call tmpmap%release(); TEST_IERR()
    call Map%release(); TEST_IERR()
    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getRangeMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getRangeMap)

  ! -------------------------------getCrsGraph-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getCrsGraph)
    type(TpetraCrsMatrix) :: mat
    type(TpetraCrsGraph) :: graph

    OUT0("Starting TpetraCrsMatrix_getCrsGraph")

    ! get identity matrix
    call TPetra_CrsMatrix_CreateIdentity(comm, mat)
    call mat%fillComplete()

    TEST_NOTHROW(graph = Mat%getCrsGraph())

    call graph%release(); TEST_IERR()
    call mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getCrsGraph")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getCrsGraph)

  ! -----------------------------getGlobalNumRows----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalNumRows)
    type(TpetraMap) :: Map
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: num_ent_per_row, ires
    integer :: lnum_local=2, commsize

    OUT0("Starting TpetraCrsMatrix_getGlobalNumRows!")

    ! create a Map
    Map = TpetraMap(TPETRA_GLOBAL_INVALID, lnum_local, comm); TEST_IERR()
    num_ent_per_row = 2
    Mat = TpetraCrsMatrix(Map, Map, num_ent_per_row, TpetraStaticProfile)
    call mat%fillComplete()

    commsize = comm%getSize()
    ires = Mat%getGlobalNumRows()
    TEST_EQUALITY(ires, int(2*commsize,size_type))

    call Map%release(); TEST_IERR()
    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getGlobalNumRows!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalNumRows)

  ! -----------------------------getGlobalNumCols----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalNumCols)
    type(TpetraMap) :: Map
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: ires
    integer :: commsize

    OUT0("Starting TpetraCrsMatrix_getGlobalNumCols!")

    commsize = comm%getSize()

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    ires = Mat%getGlobalNumCols()
    TEST_EQUALITY(ires, int(test_matrix_num_local()*commsize,size_type))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getGlobalNumCols!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalNumCols)

  ! ------------------------------getNodeNumRows------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNodeNumRows)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsMatrix_getNodeNumRows!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    ires = Mat%getNodeNumRows()
    TEST_EQUALITY(ires, int(test_matrix_num_local(),size_type))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getNodeNumRows!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNodeNumRows)

  ! ------------------------------getNodeNumCols------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNodeNumCols)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: ires
    integer :: commsize

    OUT0("Starting TpetraCrsMatrix_getNodeNumCols!")

    commsize = comm%getSize()

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    ires = Mat%getNodeNumCols()
    TEST_EQUALITY(ires, int(test_matrix_num_local(),size_type))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getNodeNumCols!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNodeNumCols)

  ! ---------------------------getGlobalNumEntries---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalNumEntries)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: ires
    integer :: commsize

    OUT0("Starting TpetraCrsMatrix_getGlobalNumEntries!")

    commsize = comm%getSize()

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    ires = Mat%getGlobalNumEntries()
    TEST_EQUALITY(ires, int(commsize*test_matrix_num_local(),size_type))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getGlobalNumEntries!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalNumEntries)

  ! ----------------------------getNodeNumEntries----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNodeNumEntries)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsMatrix_getNodeNumEntries!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    ires = Mat%getNodeNumEntries()
    TEST_EQUALITY(ires, int(test_matrix_num_local(),size_type))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getNodeNumEntries!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNodeNumEntries)

  ! -------------------------getNumEntriesInGlobalRow------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNumEntriesInGlobalRow)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: ires
    integer :: commsize, my_rank
    integer(global_ordinal_type) :: base

    OUT0("Starting TpetraCrsMatrix_getNumEntriesInGlobalRow!")

    commsize = comm%getSize()
    my_rank = comm%getRank()

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    base = test_matrix_num_local() * my_rank
    call Mat%fillComplete()

    ires = Mat%getNumEntriesInGlobalRow(int(base+1,global_ordinal_type))
    TEST_EQUALITY(ires, int(1,size_type))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getNumEntriesInGlobalRow!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNumEntriesInGlobalRow)

  ! -------------------------getNumEntriesInLocalRow-------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNumEntriesInLocalRow)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsMatrix_getNumEntriesInLocalRow!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    ires = Mat%getNumEntriesInLocalRow(1)
    TEST_EQUALITY(ires, int(1,size_type))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getNumEntriesInLocalRow!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNumEntriesInLocalRow)

  ! ------------------------getGlobalMaxNumRowEntries------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalMaxNumRowEntries)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsMatrix_getGlobalMaxNumRowEntries!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    ires = Mat%getGlobalMaxNumRowEntries()
    TEST_EQUALITY(ires, int(1,size_type))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getGlobalMaxNumRowEntries!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalMaxNumRowEntries)

  ! -------------------------getNodeMaxNumRowEntries-------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNodeMaxNumRowEntries)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: ires

    OUT0("Starting TpetraCrsMatrix_getNodeMaxNumRowEntries!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    ires = Mat%getNodeMaxNumRowEntries()
    TEST_EQUALITY(ires, int(1,size_type))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getNodeMaxNumRowEntries!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getNodeMaxNumRowEntries)

  ! -----------------------------getGlobalRowCopy----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalRowCopy)
    type(TpetraMap) :: rowMap
    integer :: lclNumRows, lclRow
    integer :: numProcs
    integer(global_ordinal_type) gblNumRows, gblNumCols, gblRow
    integer(global_ordinal_type), allocatable :: gblColInds(:)
    real(mag_type), allocatable :: vals(:)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type), parameter :: ithree=3
    integer(global_ordinal_type) :: numEnt=0

    OUT0("Starting TpetraCrsMatrix_getGlobalRowCopy!")
    call Tpetra_CrsMatrix_CreateLabelledMat(comm, Mat)
    ! Do not call this as the re-ordering done will break tests
    !call Mat%fillComplete()
    rowmap = Mat%getRowMap()

    numProcs = comm%getSize()
    gblNumRows = test_matrix_num_local() * numProcs
    gblNumCols = gblNumRows

    ! Make the arrays bigger than necessary for simplicity
    allocate(gblColInds(5),vals(5))

    do lclRow = 1, test_matrix_num_local()
      gblRow = rowMap%getGlobalElement(lclRow)
      TEST_NOTHROW(call Mat%getGlobalRowCopy(gblRow, gblColInds, vals, numEnt))

      TEST_EQUALITY(int(numEnt), 2)

      TEST_EQUALITY(gblColInds(1), 1 + modulo(gblRow + 0, gblNumCols))
      TEST_EQUALITY(gblColInds(2), 1 + modulo(gblRow + 1, gblNumCols))

      TEST_EQUALITY(int(vals(1)), int(gblColInds(1)))
      TEST_EQUALITY(int(vals(2)), int(gblColInds(2)))
    end do

    call Mat%release(); TEST_IERR()
    call rowMap%release(); TEST_IERR()
    deallocate(gblColInds,vals)

    OUT0("Finished TpetraCrsMatrix_getGlobalRowCopy!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalRowCopy)

  ! -----------------------------getLocalRowCopy------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getLocalRowCopy)
    type(TpetraMap) :: rowMap, colMap
    integer(int_type) :: lclNumRows, lclNumCols, lclRow
    integer(global_ordinal_type) gblNumRows
    integer(int_type) :: numProcs
    integer(int_type), dimension(:), allocatable :: lclColInds
    real(mag_type), allocatable :: vals(:)
    integer(int_type), dimension(:), allocatable :: curLclColInds
    real(mag_type), allocatable :: curVals(:)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type), parameter :: ithree=3
    integer(size_type) :: numEnt=0

    OUT0("Starting TpetraCrsMatrix_getLocalRowCopy!")

    ! Note that this test is currently identical to TpetraCrsMatrix_insertLocalValues
    ! because they both test getLocalRowCopy and insertLocalValues
    numProcs = comm%getSize()
    lclNumRows = 10
    gblNumRows = lclNumRows * numProcs
    rowMap = TpetraMap(gblNumRows, lclNumRows, comm)

    ! We need a locally indexed matrix to test getLocalRowCopy, so we
    ! need a column Map.  It can be the same as the row Map in this
    ! case, since we're not testing global effects.
    colMap = rowMap
    lclNumCols = lclNumRows

    ! Leave room for three locally owned entries per row.
    ! We will only fill in two entries per row.
    Mat = TpetraCrsMatrix(rowMap, colMap, ithree, TpetraStaticProfile)

    allocate(lclColInds(2))
    allocate(vals(2))

    do lclRow = 1, lclNumRows
      lclColInds(1) = 1 + modulo(lclRow + 0, lclNumCols)
      lclColInds(2) = 1 + modulo(lclRow + 1, lclNumCols)

      vals(1) = lclColInds(1)
      vals(2) = lclColInds(2)

      TEST_NOTHROW(call Mat%insertLocalValues(lclRow, lclColInds, vals))
    end do

    TEST_ASSERT(.not. Mat%isFillComplete());

    ! Make the arrays bigger than necessary, just to make sure that
    ! the methods behave correctly.
    allocate(curLclColInds(5))
    allocate(curVals(5))

    do lclRow = 1, lclNumRows
      TEST_NOTHROW(call Mat%getLocalRowCopy(lclRow, curLclColInds, curVals, numEnt))
      TEST_EQUALITY(int(numEnt), 2)

      TEST_EQUALITY(curLclColInds(1), 1 + modulo(lclRow + 0, lclNumCols));
      TEST_EQUALITY(curLclColInds(2), 1 + modulo(lclRow + 1, lclNumCols));

      TEST_EQUALITY(int(curVals(1)), curLclColInds(1));
      TEST_EQUALITY(int(curVals(2)), curLclColInds(2));
    end do

    call Mat%release(); TEST_IERR()
    call rowMap%release(); TEST_IERR()
    call colMap%release(); TEST_IERR()
    deallocate(lclColInds)
    deallocate(vals)
    deallocate(curLclColInds)
    deallocate(curVals)

    OUT0("Finished TpetraCrsMatrix_getLocalRowCopy!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getLocalRowCopy)

  ! -----------------------------getGlobalRowView----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_GetGlobalRowView)
    type(TpetraMap) :: rowMap
    type(TpetraCrsMatrix) :: Mat
    integer :: numrow, irow, nnz
    integer(global_ordinal_type), dimension(:), pointer :: view_cols
    real(scalar_type), dimension(:), pointer :: view_vals
    integer(global_ordinal_type) :: gblrow
    integer(global_ordinal_type), allocatable :: compare_cols(:)
    real(scalar_type), allocatable :: compare_vals(:)

    OUT0("Starting TpetraCrsMatrix_GetGlobalRowView")

    call TPetra_CrsMatrix_CreateTestMatrix_A(comm,Mat)

    numrow = Mat%getNodeNumRows()
    rowmap = Mat%getRowMap()
    do irow=1,numrow
       gblrow = rowmap%getGlobalElement(irow)
       TEST_NOTHROW(call Mat%getGlobalRowView(gblrow, view_cols, view_vals))

       ! Check the results row by row
       call Tpetra_CrsMatrix_GetTestMatrixRow_A(comm, irow, gblrow, compare_cols, compare_vals, nnz)
       TEST_FLOATING_ARRAY_EQUALITY(view_vals, compare_vals, epsilon(0.d0))
       TEST_ARRAY_EQUALITY(view_cols, compare_cols)
       deallocate(compare_cols, compare_vals)
    enddo

    nullify(view_cols, view_vals)
    call Mat%release();  TEST_IERR()
    call rowmap%release();  TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_GetGlobalRowView")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_GetGlobalRowView)

  ! ------------------------------getProfileType------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getProfileType)
    type(TpetraCrsMatrix) :: Mat

    OUT0("Starting TpetraCrsMatrix_getProfileType!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)

    ! This test fails if fillComplete is called
    !call Mat%fillComplete()

    TEST_ASSERT(Mat%getProfileType() == TpetraStaticProfile)

    call Mat%release();  TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getProfileType!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getProfileType)

  ! -----------------------------getFrobeniusNorm----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getFrobeniusNorm)
    type(TpetraCrsMatrix) :: Mat
    real(mag_type) :: fresult, fnorm

    OUT0("Starting TpetraCrsMatrix_getFrobeniusNorm!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)

    fnorm = sqrt(real(comm%getSize()*test_matrix_num_local(), kind=scalar_type))
    TEST_NOTHROW(fresult=Mat%getFrobeniusNorm())
    TEST_ASSERT(fresult==fnorm)

    call Mat%release();  TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getFrobeniusNorm!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getFrobeniusNorm)

  ! -------------------------------getAllValues------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getAllValues)
    type(TpetraMap) :: colmap
    type(TpetraCrsMatrix) :: Mat
    integer :: numrow, numnnz, ii
    integer(size_type), dimension(:), allocatable :: row_ptrs, rp_res
    integer(int_type), dimension(:), allocatable :: colind, col_res, col_gbl
    real(scalar_type), dimension(:), allocatable :: values, val_res

    OUT0("Starting TpetraCrsMatrix_getAllValues!")

    call TPetra_CrsMatrix_CreateTestMatrix_A(comm, Mat)
    call Mat%fillComplete()

    numrow = Mat%getNodeNumRows()
    numnnz = Mat%getNodeNumEntries()
    allocate(row_ptrs(numrow+1),rp_res(numrow+1))
    allocate(colind(numnnz),col_res(numnnz),col_gbl(numnnz))
    allocate(values(numnnz),val_res(numnnz))

    call Mat%getAllValues(row_ptrs, colind, values)

    ! getAllValues returns the local indexing of column -- convert to global
    ! for the comparison
    colmap=Mat%getColMap()
    do ii = 1,numnnz
       col_gbl(ii) = colmap%getGlobalElement(colind(ii))
    enddo

    !---------------------------------------------------------------
    ! CHECK results
    !---------------------------------------------------------------
    ! Get expected results
    call TPetra_CrsMatrix_getTestMatrixResults(comm, rp_res, col_res, val_res)

    TEST_ARRAY_EQUALITY(row_ptrs, rp_res)
    TEST_FLOATING_ARRAY_EQUALITY(values, val_res, epsilon(0.d0))

    ! column indices are often reordered so test the sum's of the global indices
    ! this is why global indices are preferred -- larger sums less likely if err
    TEST_EQUALITY(SUM(col_gbl), SUM(col_res))

    deallocate(row_ptrs, colind, values, rp_res, col_res, val_res, col_gbl)

    call colmap%release();  TEST_IERR()
    call Mat%release();  TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_getAllValues!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getAllValues)

  ! -------------------------------setAllValues------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_setAllValues)
    type(TpetraCrsMatrix) :: Mat, mat_compare
    type(TpetraMap) :: map, colmap
    integer :: numrow, numnnz, ii, id
    integer(global_ordinal_type), dimension(:), allocatable :: rowptr, rowchk
    integer, dimension(:), allocatable :: ind, indchk, indgbl, indlcl
    real(mag_type), dimension(:), allocatable :: vals, valchk

    OUT0("Starting TpetraCrsMatrix_setAllValues!")

    map = TpetraMap(TPETRA_GLOBAL_INVALID, test_matrix_num_row(), comm)
    Mat = TpetraCrsMatrix(map, map, test_matrix_max_entries_per_row(), TpetraStaticProfile)

    ! Like getAllValues use the TestMatrix.
    call TPetra_CrsMatrix_CreateTestMatrix_A(comm, mat_compare)
    call mat_compare%fillComplete()

    numrow = mat_compare%getNodeNumRows()
    numnnz = mat_compare%getNodeNumEntries()
    allocate(rowptr(numrow+1), rowchk(numrow+1))
    allocate(indgbl(numnnz), ind(numnnz), indchk(numnnz), indlcl(numnnz))
    allocate(vals(numnnz), valchk(numnnz))

    call TPetra_CrsMatrix_getTestMatrixResults(comm, rowchk, indchk, valchk)

    ! setAllValues wants local indicies; helper returns globals so use ColMap to convert
    colmap=Mat%getColMap()
    do ii = 1,numnnz
       indlcl(ii) = colmap%getLocalElement(int(indchk(ii), global_ordinal_type))
    enddo

    TEST_NOTHROW(call Mat%setAllValues(rowchk, indlcl, valchk))


    !---------------------------------------------------------------
    ! CHECK results
    !---------------------------------------------------------------
    call Mat%getAllValues(rowptr, ind, vals)
    TEST_ARRAY_EQUALITY(rowptr, rowchk)
    TEST_FLOATING_ARRAY_EQUALITY(vals, valchk, epsilon(0.d0))

    deallocate(rowptr, ind, vals, rowchk, indchk, valchk, indgbl, indlcl)
    call Mat%release(); TEST_IERR()
    call colmap%release(); TEST_IERR()
    call Map%release(); TEST_IERR()
    call mat_compare%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_setAllValues!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_setAllValues)

  ! ----------------------------insertGlobalValues---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_insertGlobalValues)
    type(TpetraMap) :: colmap
    type(TpetraCrsMatrix) :: Mat
    integer :: numnnz, ii
    integer :: numrow = 4
    integer(size_type), dimension(:), allocatable :: row_ptrs, rp_res
    integer(int_type), dimension(:), allocatable :: colind, col_res, col_gbl
    real(scalar_type), dimension(:), allocatable :: values, val_res
    type(TpetraMap) :: map
    integer :: nnz, irow
    integer(global_ordinal_type) :: gblrow
    integer(global_ordinal_type), allocatable :: cols(:)
    real(scalar_type), allocatable :: vals(:)

    OUT0("Starting TpetraCrsMatrix_insertGlobalValues!")

    ! create a Map and matrix
    map = TpetraMap(TPETRA_GLOBAL_INVALID, numrow, comm)

    Mat = TpetraCrsMatrix(map, test_matrix_max_entries_per_row(), TpetraStaticProfile)

    ! The test setup TPetra_CrsMatrix_CreateTestMatrix is essentially
    ! a test of insertGlobalValues but we'd like to have it appear
    ! explicitly in this test so we duplicate the code here.
    do irow=1,numrow
      call Tpetra_CrsMatrix_GetTestMatrixRow_A(comm, irow, gblrow, cols, vals, nnz)
      call Mat%insertGlobalValues(gblrow, cols, vals)
      deallocate(cols,vals)
    enddo

    !---------------------------------------------------------------
    ! Now read values and check them
    !---------------------------------------------------------------
    call Mat%fillComplete()

    numrow = Mat%getNodeNumRows()
    numnnz = Mat%getNodeNumEntries()
    allocate(row_ptrs(numrow+1),rp_res(numrow+1))
    allocate(colind(numnnz),col_res(numnnz),col_gbl(numnnz))
    allocate(values(numnnz),val_res(numnnz))

    call Mat%getAllValues(row_ptrs, colind, values)

    ! getAllValues returns the local indexing of column -- convert to global
    ! for the comparison
    colmap=Mat%getColMap()
    do ii = 1,numnnz
       col_gbl(ii) = colmap%getGlobalElement(colind(ii))
    enddo

    !---------------------------------------------------------------
    ! CHECK results
    !---------------------------------------------------------------
    ! Get expected results
    call TPetra_CrsMatrix_getTestMatrixResults(comm, rp_res, col_res, val_res)

    TEST_ARRAY_EQUALITY(row_ptrs, rp_res)
    TEST_FLOATING_ARRAY_EQUALITY(values, val_res, epsilon(0.d0))

    ! column indices are often reordered so test the sum's of the global indices
    ! this is why global indices are preferred -- larger sums less likely if err
    TEST_EQUALITY(SUM(col_gbl), SUM(col_res))

    deallocate(row_ptrs, colind, values, rp_res, col_res, val_res, col_gbl)
    call map%release();  TEST_IERR()
    call colmap%release();  TEST_IERR()
    call Mat%release();  TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_insertGlobalValues!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_insertGlobalValues)

  ! ----------------------------insertLocalValues----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_insertLocalValues)
    type(TpetraMap) :: rowMap, colMap
    integer(int_type) :: lclNumRows, lclNumCols, lclRow
    integer(global_ordinal_type) gblNumRows
    integer(int_type) :: numProcs
    integer(int_type), dimension(:), allocatable :: lclColInds
    real(mag_type), allocatable :: vals(:)
    integer(int_type), dimension(:), allocatable :: curLclColInds
    real(mag_type), allocatable :: curVals(:)
    type(TpetraCrsMatrix) :: Mat
    integer(size_type) :: numEnt=0

    ! Note that this test is currently identical to TpetraCrsMatrix_getLocalRowCopy
    ! becuase they both test getLocalRowCopy and insertLocalValues

    OUT0("Starting TpetraCrsMatrix_insertLocalValues!")

    numProcs = comm%getSize()
    lclNumRows = 10
    gblNumRows = lclNumRows * numProcs
    rowMap = TpetraMap(gblNumRows, lclNumRows, comm)

    ! We need a locally indexed matrix to test getLocalRowCopy, so we
    ! need a column Map.  It can be the same as the row Map in this
    ! case, since we're not testing global effects.
    colMap = rowMap
    lclNumCols = lclNumRows

    ! Leave room for three locally owned entries per row.
    ! We will only fill in two entries per row.
    Mat = TpetraCrsMatrix(rowMap, colMap, test_matrix_max_entries_per_row(), TpetraStaticProfile)

    allocate(lclColInds(2))
    allocate(vals(2))

    do lclRow = 1, lclNumRows
      lclColInds(1) = 1 + modulo(lclRow + 0, lclNumCols)
      lclColInds(2) = 1 + modulo(lclRow + 1, lclNumCols)

      vals(1) = lclColInds(1)
      vals(2) = lclColInds(2)

      TEST_NOTHROW(call Mat%insertLocalValues(lclRow, lclColInds, vals))
    end do

    TEST_ASSERT(.not. Mat%isFillComplete());

    ! Make the arrays bigger than necessary, just to make sure that
    ! the methods behave correctly.
    allocate(curLclColInds(5))
    allocate(curVals(5))

    do lclRow = 1, lclNumRows
      TEST_NOTHROW(call Mat%getLocalRowCopy(lclRow, curLclColInds, curVals, numEnt))
      TEST_EQUALITY(int(numEnt), 2)

      TEST_EQUALITY(curLclColInds(1), 1 + modulo(lclRow + 0, lclNumCols));
      TEST_EQUALITY(curLclColInds(2), 1 + modulo(lclRow + 1, lclNumCols));

      TEST_EQUALITY(int(curVals(1)), curLclColInds(1));
      TEST_EQUALITY(int(curVals(2)), curLclColInds(2));
    end do

    call Mat%release(); TEST_IERR()
    call rowMap%release(); TEST_IERR()
    call colMap%release(); TEST_IERR()
    deallocate(lclColInds)
    deallocate(vals)
    deallocate(curLclColInds)
    deallocate(curVals)

    OUT0("Finished TpetraCrsMatrix_insertLocalValues!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_insertLocalValues)

  ! ---------------------------replaceGlobalValues---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceGlobalValues)
    type(TpetraCrsMatrix) :: Mat
    integer :: irow, numvalid
    integer(global_ordinal_type) :: base, gblrow, cols(1)
    real(mag_type) :: vals(1)=[2.]
    integer :: lclrow

    OUT0("Starting TpetraCrsMatrix_replaceGlobalValues!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)

    base = test_matrix_num_local() * comm%getRank()
    do irow = 1, test_matrix_num_local()
      gblrow = base + int(irow, kind=global_ordinal_type)
      cols(1) = gblrow
      TEST_NOTHROW(numvalid = Mat%replaceGlobalValues(gblrow, cols, vals))
      ! Only one entry changes per iteration, and numvalid indicates how many
      ! substitutions have occured
      TEST_EQUALITY(numvalid, 1)
    end do

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_replaceGlobalValues!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceGlobalValues)

  ! ------------------------------setAllToScalar------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_setAllToScalar)
    type(TpetraCrsMatrix) :: Mat
    type(TpetraMap) :: rowmap
    integer :: lclRow, ii
    real(mag_type), allocatable :: vals(:)
    integer(global_ordinal_type) :: gblRow, numEnt = 0
    integer(global_ordinal_type), allocatable :: gblColInds(:)

    OUT0("Starting TpetraCrsMatrix_setAllToScalar!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    rowmap = Mat%getRowMap()

    TEST_NOTHROW(call Mat%setAllToScalar(real(2,mag_type)))
    call Mat%fillComplete()
    allocate(gblColInds(test_matrix_num_local()), vals(test_matrix_num_local()))

    do lclRow = 1, test_matrix_num_local()
      gblRow = rowMap%getGlobalElement(lclRow)
      TEST_NOTHROW(call Mat%getGlobalRowCopy(gblRow, gblColInds, vals, numEnt))

      do ii = 1, numEnt
        TEST_FLOATING_EQUALITY(vals(ii), real(2, mag_type), epsilon(real(2, mag_type)))
      end do
    end do

    call Mat%release(); TEST_IERR()
    call rowmap%release(); TEST_IERR()
    deallocate(vals, gblColInds)
    OUT0("Finished TpetraCrsMatrix_setAllToScalar!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_setAllToScalar)

  ! ----------------------------------scale----------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_scale)
    type(TpetraCrsMatrix) :: Mat
    real(mag_type) :: norm, norm_s

    OUT0("Starting TpetraCrsMatrix_scale!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()
    norm = Mat%getFrobeniusNorm()

    call Mat%resumeFill()
    TEST_NOTHROW(call Mat%scale(real(2.0, mag_type)))
    call Mat%fillComplete()
    norm_s = Mat%getFrobeniusNorm()

    !Norm properties: Scaling matrix by 2 should scale norm by 2
    TEST_FLOATING_EQUALITY(norm_s, 2.0_mag_type * norm, epsilon(norm))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_scale!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_scale)

  ! ------------------------------globalAssemble------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_globalAssemble)
    type(TpetraCrsMatrix) :: Mat

    OUT0("Starting TpetraCrsMatrix_globalAssemble!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    TEST_THROW(call Mat%globalAssemble())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_globalAssemble!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_globalAssemble)

  ! ------------------------------replaceColMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceColMap)
    type(TpetraMap) :: map, newcolmap
    type(TpetraCrsMatrix) :: mat
    integer(size_type), parameter :: ione=1

    OUT0("Starting TpetraCrsMatrix_replaceColMap")

    map = TpetraMap(TPETRA_GLOBAL_INVALID, test_matrix_num_local(), comm); TEST_IERR()
    mat = TpetraCrsMatrix(map, ione, TpetraStaticProfile)

    newcolmap = TpetraMap(TPETRA_GLOBAL_INVALID, test_matrix_num_local()-1, comm); TEST_IERR()

    call mat%replaceColMap(newcolmap); TEST_IERR()

    call mat%release(); TEST_IERR()
    call map%release(); TEST_IERR()
    call newcolmap%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_replaceColMap")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceColMap)

  ! -----------------------removeEmptyProcessesInPlace------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_removeEmptyProcessesInPlace)
    type(TpetraCrsMatrix) :: mat
    type(TpetraMap) :: map
    integer(size_type), parameter :: ione=1
    integer :: lnum_local=3

    OUT0("Starting TpetraCrsMatrix_removeEmptyProcessesInPlace")

    map = TpetraMap(TPETRA_GLOBAL_INVALID, lnum_local, comm); TEST_IERR()
    mat = TpetraCrsMatrix(map, ione, TpetraStaticProfile)

    call mat%removeEmptyProcessesInPlace(map); TEST_IERR()

    call map%release(); TEST_IERR()
    call mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_removeEmptyProcessesInPlace")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_removeEmptyProcessesInPlace)

  ! --------------------------------hasColMap--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_hasColMap)
    type(TpetraCrsMatrix) :: Mat

    OUT0("Starting TpetraCrsMatrix_hasColMap!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    TEST_ASSERT(Mat%hasColMap())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_hasColMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_hasColMap)

  ! -----------------------------isLocallyIndexed----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isLocallyIndexed)
    type(TpetraCrsMatrix) :: Mat

    OUT0("Starting TpetraCrsMatrix_isLocallyIndexed!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    TEST_ASSERT(Mat%isLocallyIndexed())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_isLocallyIndexed!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isLocallyIndexed)

  ! ----------------------------isGloballyIndexed----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isGloballyIndexed)
    type(TpetraCrsMatrix) :: Mat

    OUT0("Starting TpetraCrsMatrix_isGloballyIndexed!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)

    TEST_ASSERT(Mat%isGloballyIndexed())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_isGloballyIndexed!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isGloballyIndexed)

  ! ------------------------------isFillComplete------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isFillComplete)
    type(TpetraCrsMatrix) :: Mat

    OUT0("Starting TpetraCrsMatrix_isFillComplete!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    TEST_ASSERT(Mat%isFillComplete())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_isFillComplete!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isFillComplete)

  ! -------------------------------isFillActive------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isFillActive)
    type(TpetraCrsMatrix) :: Mat

    OUT0("Starting TpetraCrsMatrix_isFillActive!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)

    TEST_ASSERT(Mat%isFillActive())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_isFillActive!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isFillActive)

  ! ----------------------------isStorageOptimized---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isStorageOptimized)
    type(TpetraCrsMatrix) :: mat

    OUT0("Starting TpetraCrsMatrix_isStorageOptimized")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    TEST_ASSERT(mat%isStorageOptimized())

    call mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_isStorageOptimized")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isStorageOptimized)

  ! ------------------------------isStaticGraph------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isStaticGraph)
    type(TpetraMap) :: map
    type(TpetraCrsMatrix) :: mat
    type(TpetraCrsGraph) :: graph
    integer(size_type), parameter :: ione=1

    OUT0("Starting TpetraCrsMatrix_isStaticGraph")

    map = TpetraMap(TPETRA_GLOBAL_INVALID, test_matrix_num_local(), comm); TEST_IERR()
    graph = TpetraCrsGraph(Map, Map, ione, TpetraStaticProfile)
    call graph%fillComplete()
    mat = TpetraCrsMatrix(graph)

    TEST_ASSERT(mat%isStaticGraph())

    call map%release(); TEST_IERR()
    call mat%release(); TEST_IERR()
    call graph%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_isStaticGraph")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_isStaticGraph)

  ! -----------------------------supportsRowViews----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_supportsRowViews)
    type(TpetraCrsMatrix) :: Mat

    OUT0("Starting TpetraCrsMatrix_supportsRowViews")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    TEST_ASSERT(Mat%supportsRowViews())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_supportsRowViews")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_supportsRowViews)

  ! ----------------------------hasTransposeApply----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_hasTransposeApply)
    type(TpetraCrsMatrix) :: Mat

    OUT0("Starting TpetraCrsMatrix_hasTransposeApply!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    TEST_ASSERT(Mat%hasTransposeApply())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_hasTransposeApply!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_hasTransposeApply)

  ! -------------------------------description-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_description)
    type(TpetraCrsMatrix) :: Mat
    character(kind=C_CHAR, len=:), allocatable :: description
    character(kind=C_CHAR, len=200) :: desc_res
    character(kind=C_CHAR, len=2) :: str_size
    integer :: num_images

    OUT0("Starting TpetraCrsMatrix_description!")

    num_images = comm%getSize()
    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    TEST_NOTHROW(description = Mat%description())

    ! Test against known result
    ! string processing in fortran:
    !  FORTRAN was the language of choice for the same reason that three-legged races are popular. - -  Kernighan
    write(str_size,fmt='(i2.2)') num_images*10
    desc_res="Tpetra::CrsMatrix (Kokkos refactor): {isFillComplete: true, global dimensions: ["
    desc_res=TRIM(desc_res)//str_size//', '//str_size//'], global number of entries: '//str_size//'}'
    TEST_EQUALITY(TRIM(description),TRIM(desc_res))

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_description!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_description)

  ! --------------------------computeGlobalConstants-------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_computeGlobalConstants)
    type(TpetraCrsMatrix) :: Mat
    character(kind=C_CHAR, len=:), allocatable :: description

    OUT0("Starting TpetraCrsMatrix_computeGlobalConstants!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    TEST_NOTHROW(call Mat%computeGlobalConstants())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_computeGlobalConstants!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_computeGlobalConstants)

  ! ---------------------------haveGlobalConstants---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_haveGlobalConstants)
    type(TpetraCrsMatrix) :: Mat
    character(kind=C_CHAR, len=:), allocatable :: description

    OUT0("Starting TpetraCrsMatrix_haveGlobalConstants!")

    call TPetra_CrsMatrix_CreateIdentity(comm, Mat)
    call Mat%fillComplete()

    call Mat%computeGlobalConstants()
    TEST_ASSERT(Mat%haveGlobalConstants())

    call Mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_haveGlobalConstants!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_haveGlobalConstants)

  ! -----------------------replaceDomainMapAndImporter------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceDomainMapAndImporter)
    type(TpetraCrsMatrix) :: mat
    type(TpetraMap) :: newdomainmap, newcolmap
    type(TpetraImport) :: newimporter
    integer(global_ordinal_type) :: num_global

    OUT0("Starting TpetraCrsMatrix_replaceDomainMapAndImporter")

    call TPetra_CrsMatrix_CreateIdentity(comm, mat)
    call mat%fillComplete()

    num_global = 5*comm%getSize()
    newdomainmap = TpetraMap(num_global, comm); TEST_IERR()
    newcolmap =  mat%getColMap(); TEST_IERR()
    newimporter = TpetraImport(newdomainmap, newcolmap ); TEST_IERR()

    TEST_NOTHROW(call mat%replaceDomainMapAndImporter(newdomainmap, newimporter))

    call newimporter%release(); TEST_IERR()
    call newdomainmap%release(); TEST_IERR()
    call newcolmap%release(); TEST_IERR()
    call mat%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_replaceDomainMapAndImporter")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceDomainMapAndImporter)

  ! -------------------------------gaussSeidel-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_gaussSeidel)
    type(TpetraMap) :: Map, row_map, tmp_map
    type(TeuchosComm) :: tcomm
    type(TpetraCrsMatrix) :: Mat
    type(TpetraMultiVector) :: b
    type(TpetraMultiVector) :: x
    type(TpetraMultiVector) :: d
    real(scalar_type) :: dampingfactor
    real(scalar_type), parameter :: one=1.
    integer(size_type), parameter :: ione=1
    integer(size_type), parameter :: num_vecs=1
    real(mag_type) :: vals(1)=[1.], normb(num_vecs), normx(num_vecs)
    integer(size_type) :: num_images, my_image_id
    integer :: irow
    integer(global_ordinal_type) :: base, gblrow, cols(1)
    integer(kind(TpetraESweepDirection)) :: direction
    integer(C_INT) :: numsweeps

    OUT0("Starting TpetraCrsMatrix_gaussSeidel")

    my_image_id = comm%getRank()

    ! create a Map
    Map = TpetraMap(TPETRA_GLOBAL_INVALID, test_matrix_num_local(), comm); TEST_IERR()

    ! create the identity matrix
    base = test_matrix_num_local() * my_image_id
    Mat = TpetraCrsMatrix(Map, ione, TpetraStaticProfile)
    do irow = 1, test_matrix_num_local()
      gblrow = base + int(irow, kind=global_ordinal_type)
      cols(1) = gblrow
      vals(1) = one
      call Mat%insertGlobalValues(gblrow, cols, vals)
    end do
    call mat%fillComplete()

    b = TpetraMultiVector(Map, num_vecs); TEST_IERR()
    call b%randomize(); TEST_IERR()
    x = TpetraMultiVector(Map, num_vecs); TEST_IERR()
    call x%randomize(); TEST_IERR()
    d = TpetraMultiVector(Map, num_vecs); TEST_IERR()
    call d%putScalar(-1.0_scalar_type)  ! Inverse diagonal
    dampingfactor = 0.04
    direction = 0
    numsweeps = 4

    call b%norm2(normb)
    call x%norm2(normx)
    !write(*,*) 'Before iterating: normx= ', normx,' normb ', normb

    call Mat%gaussSeidel(b, x, d, dampingfactor, direction, numsweeps); TEST_IERR()

    call b%norm2(normb)
    call x%norm2(normx)
    !write(*,*) 'After iterating: normx= ', normx,' normb ', normb

    call b%release(); TEST_IERR()
    call x%release(); TEST_IERR()
    call d%release(); TEST_IERR()
    call Mat%release(); TEST_IERR()
    call Map%release(); TEST_IERR()

    OUT0("Finished TpetraCrsMatrix_gaussSeidel")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_gaussSeidel)

  ! -----------------------------gaussSeidelCopy------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_gaussSeidelCopy)
    type(TpetraCrsMatrix) :: Mat
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
    !Mat = TpetraCrsMatrix(); TEST_IERR()
    !call Mat%gaussSeidelCopy(x, b, d, dampingfactor, direction, numsweeps, zeroinitialguess); TEST_IERR()

    !call x%release(); TEST_IERR()
    !call b%release(); TEST_IERR()
    !call d%release(); TEST_IERR()
    !call Mat%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_gaussSeidelCopy: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_gaussSeidelCopy")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_gaussSeidelCopy)

end program test_tpetraCrsMatrix
