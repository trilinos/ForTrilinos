program test_TpetraCrsMatrix
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
#endif

  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_Basic1)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_AlphaBetaMultiply)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_ActiveFill)
  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_SimpleEigTest)

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
!  ADD_SUBTEST_AND_RUN(TpetraCrsMatrix_getGlobalRowCopy)

  call comm%release()
  CHECK_IERR()

  SHUTDOWN_TEST()

contains

  ! ------------------------------- Basic1 ----------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_Basic1)
    type(TpetraMap) :: Map, row_map
    type(TeuchosComm) :: tcomm
    type(TpetraCrsMatrix) :: Mat
    type(TpetraMultiVector) :: mvrand, mvres
    type(string) :: description
    logical(c_bool), parameter :: false=.false., true=.true.
    integer(size_type) :: num_images, my_image_id
    integer(size_type), parameter :: num_local=10, num_vecs=5
    integer(size_type), parameter :: ione=1
    integer(local_ordinal_type) :: irow
    integer(global_ordinal_type) :: base, gblrow, cols(1)
    real(mag_type) :: vals(1)=[1.], norms(num_vecs), zeros(num_vecs), fnorm
    real(scalar_type), parameter :: zero=0., one=1., negone=-1., two=2., four=4.

    type(TeuchosArrayViewDouble) :: av1, av2
    real(scalar_type), allocatable :: a1(:), a2(:)
    integer(size_type) :: lda

    OUT0("Starting TpetraCrsMatrix_Basic1")

    zeros = zero

    num_images = comm%getSize()
    my_image_id = comm%getRank()

    ! create a Map
    call Map%create(invalid, num_local, index_base, comm); TEST_IERR()
    call mvrand%create(Map, num_vecs, false); TEST_IERR()
    call mvres%create(Map, num_vecs, false); TEST_IERR()
    call mvrand%randomize(); TEST_IERR()

    ! create the identity matrix
    base = num_local * my_image_id;
    call Mat%create(Map, ione, DynamicProfile)
    do irow = 1, num_local
      gblrow = base + int(irow, kind=global_ordinal_type)
      cols(1) = gblrow
      vals(1) = one
      call Mat%insertGlobalValues(gblrow, cols, vals)
    end do

    TEST_ASSERT(Mat%isGloballyIndexed())
    TEST_ASSERT((.not. Mat%isLocallyIndexed()))
    TEST_ASSERT(Mat%getProfileType() == DynamicProfile)
    call Mat%fillComplete(); TEST_IERR()
    row_map = Mat%getRowMap()

    ! test the properties
    TEST_EQUALITY(Mat%getGlobalNumEntries(), num_images*num_local)
    TEST_EQUALITY(Mat%getNodeNumEntries(), num_local)
    TEST_EQUALITY(Mat%getGlobalNumRows(), num_images*num_local)
    TEST_EQUALITY(Mat%getGlobalNumCols(), num_images*num_local)
    TEST_EQUALITY(Mat%getNodeNumRows(), num_local)
    TEST_EQUALITY(Mat%getNodeNumCols(), num_local)
    TEST_EQUALITY(Mat%getGlobalNumDiags(), num_images*num_local)
    TEST_EQUALITY(Mat%getNodeNumDiags(), num_local)
    TEST_EQUALITY(Mat%getGlobalMaxNumRowEntries(), 1)
    TEST_EQUALITY(Mat%getNodeMaxNumRowEntries(), 1)
    TEST_ASSERT(Mat%isFillComplete())
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_EQUALITY(Mat%getIndexBase(), 1)
    fnorm = sqrt(real(num_images*num_local, kind=scalar_type))
    TEST_EQUALITY(Mat%getFrobeniusNorm(), fnorm)
    TEST_ASSERT(row_map%isSameAs(Mat%getColMap()))
    TEST_ASSERT(row_map%isSameAs(Mat%getDomainMap()))
    TEST_ASSERT(row_map%isSameAs(Mat%getRangeMap()))
    TEST_ASSERT(Mat%hasColMap())
    TEST_ASSERT(Mat%haveGlobalConstants())
    !TEST_ASSERT((.not. Mat%isLowerTriangular())) TODO: This throws an error from Tpetra
    !TEST_ASSERT((.not. Mat%isUpperTriangular())) TODO: This tthrows an error from Tpetra
    TEST_ASSERT(Mat%isLocallyIndexed())
    TEST_ASSERT(Mat%hasTransposeApply())
    TEST_ASSERT((.not. Mat%isGloballyIndexed()))

    tcomm = Mat%getComm(); TEST_IERR()
    TEST_EQUALITY(tcomm%getRank(), comm%getRank())
    TEST_EQUALITY(tcomm%getSize(), comm%getSize())

    ! get the description, just to see if it does not throw
    description = Mat%description(); TEST_IERR()

    do irow = 1, num_local
      gblrow = row_map%getGlobalElement(irow)
      TEST_ASSERT(gblrow == (base + int(irow, kind=global_ordinal_type)))
      TEST_EQUALITY(Mat%getNumEntriesInGlobalRow(gblrow), 1)
    end do

    ! test the action
    call mvres%randomize(); TEST_IERR()
    call Mat%apply(mvrand, mvres); TEST_IERR()
    call mvres%update(one, mvrand, negone); TEST_IERR()

    call mvres%norm1(norms); TEST_IERR()
    TEST_COMPARE_FLOATING_ARRAYS(norms, zero, epsilon(zero))

    ! Set all diagonal entries to 2 and do again
    call Mat%resumeFill(); TEST_IERR()
    call Mat%setAllToScalar(two); TEST_IERR()
    call Mat%fillComplete(); TEST_IERR()

    call mvres%randomize(); TEST_IERR()
    call Mat%apply(mvrand, mvres); TEST_IERR()
    call mvres%update(two, mvrand, negone); TEST_IERR()

    call mvres%norm1(norms); TEST_IERR()
    TEST_COMPARE_FLOATING_ARRAYS(norms, zero, epsilon(zero))

    ! Scale diagonal entries by 2 and do again
    call Mat%resumeFill(); TEST_IERR()
    call Mat%scale(two); TEST_IERR()
    call Mat%fillComplete(); TEST_IERR()

    call mvres%randomize(); TEST_IERR()
    call Mat%apply(mvrand, mvres); TEST_IERR()
    call mvres%update(four, mvrand, negone); TEST_IERR()

    call mvres%norm1(norms); TEST_IERR()
    TEST_COMPARE_FLOATING_ARRAYS(norms, zero, epsilon(zero))

    call Mat%release()
    call Map%release()
    call mvres%release()
    call mvrand%release()

    OUT0("Finished TpetraCrsMatrix_Basic1!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_Basic1)

  ! -------------------------- SimpleEigTest --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_SimpleEigTest)
    type(TpetraMap) :: Map, row_map
    type(TpetraCrsMatrix) :: A
    type(TpetraMultiVector) :: ones, threes
    integer(size_type), parameter :: izero=0, ione=1, ithree=3
    integer(size_type) :: num_images, my_image_id
    integer(local_ordinal_type) :: nnz
    logical(c_bool), parameter :: false=.false., true=.true.
    real(scalar_type), parameter :: zero=0., one=1., two=2., negthree=-3.
    real(scalar_type) :: vals(3)
    real(norm_type) :: norms(1)
    integer(global_ordinal_type) :: gblrow, cols(3)

    OUT0("Starting TpetraCrsMatrix_SimpleEigTest")
    num_images = comm%getSize()
    my_image_id = comm%getRank()

    if (num_images < 2) return

    ! create a Map
    call Map%create(invalid, ione, index_base, comm); TEST_IERR()

    ! create a multivector ones(n,1)
    call ones%create(map, ione, false); TEST_IERR()
    call threes%create(map, ione, false); TEST_IERR()
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

    call A%create(map, izero, DynamicProfile); TEST_IERR()
    gblrow = my_image_id + 1
    if (gblrow == 1) then
      nnz = 2
      cols(1:nnz) = [gblrow, gblrow+1]
      vals(1:nnz) = [two, one]
      call A%insertGlobalValues(gblrow, nnz, vals, cols); TEST_IERR()
    else if (gblrow == num_images) then
      nnz = 2;
      cols(1:nnz) = [gblrow-1, gblrow]
      vals(1:nnz) = [one, two]
      call A%insertGlobalValues(gblrow, nnz, vals, cols); TEST_IERR()
    else
      nnz = 3;
      vals = [one, one, one]
      cols = [gblrow-1, gblrow, gblrow+1]
      call A%insertGlobalValues(gblrow, nnz, vals, cols); TEST_IERR()
    end if

    call A%fillComplete(); TEST_IERR()
    row_map = A%getRowMap(); TEST_IERR()

    ! test the properties
    TEST_EQUALITY(A%getGlobalNumEntries(), 3*num_images-2)
    TEST_EQUALITY(A%getNodeNumEntries(), nnz)
    TEST_EQUALITY(A%getGlobalNumRows(), num_images)
    TEST_EQUALITY(A%getNodeNumRows(), 1)
    TEST_EQUALITY(A%getNodeNumCols(), nnz)
    TEST_EQUALITY(A%getGlobalNumDiags(), num_images)
    TEST_EQUALITY(A%getNodeNumDiags(), 1)
    if (num_images > 2) then
      TEST_EQUALITY(A%getGlobalMaxNumRowEntries(), 3)
    else
      TEST_EQUALITY(A%getGlobalMaxNumRowEntries(), 2)
    end if
    TEST_EQUALITY(A%getNodeMaxNumRowEntries(), nnz)
    TEST_EQUALITY(A%getIndexBase(), 1)
    TEST_ASSERT((.not. row_map%isSameAs(A%getColMap())))
    TEST_ASSERT(row_map%isSameAs(A%getDomainMap()))
    TEST_ASSERT(row_map%isSameAs(A%getRangeMap()))

    ! test the action
    call threes%randomize(); TEST_IERR()
    call A%apply(ones, threes); TEST_IERR()

    ! now, threes should be 3*ones
    call threes%update(negthree, ones, one)
    !call threes%norm1(norms) TODO: This call hangs with numprocs > 1
    TEST_COMPARE_FLOATING_ARRAYS(norms, zero, epsilon(zero))

    call ones%release()
    call threes%release()
    call Map%release()
    call A%release()
    call row_map%release()

    OUT0("Finished TpetraCrsMatrix_SimpleEigTest!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_SimpleEigTest)

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
    call Map%create(invalid, ithree, index_base, comm); TEST_IERR()

    ! Create the identity matrix, three rows per proc
    base = 3 * my_image_id;
    call Mat%create(Map, ione, DynamicProfile); TEST_IERR()
    do i = 1, 3
      gblrow = base + i
      cols(1) = gblrow
      call Mat%insertGlobalValues(gblrow, cols, ones); TEST_IERR()
    end do
    call Mat%fillComplete(); TEST_IERR()

    call X%create(Map, numvecs); TEST_IERR()
    call Y%create(Map, numvecs); TEST_IERR()
    call Z%create(map, numvecs); TEST_IERR()

    call random_number(alpha)
    call random_number(beta)

    call X%randomize()
    call Y%randomize()

    ! Z = alpha*X + beta*Y
    call Z%update(alpha, X, beta, Y, zero)

    ! test the action: Y = alpha*I*X + beta*Y = alpha*X + beta*Y = Z
    call Mat%apply(X, Y, NO_TRANS, alpha, beta)
    !
    call Z%norm1(normz)
    call Y%norm1(normy)
    TEST_COMPARE_FLOATING_ARRAYS(normy, normz, epsilon(zero))

    call Z%release()
    call Y%release()
    call X%release()
    call Mat%release()
    call Map%release()

    OUT0("Finished TpetraCrsMatrix_AlphaBetaMultiply!")
  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_AlphaBetaMultiply)

  ! ----------------------------- ActiveFill --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_ActiveFill)
    type(TpetraMap) :: Map
    type(ParameterList) :: params
    type(TpetraCrsMatrix) :: Mat
    integer(size_type), parameter :: izero=0, ione=1
    logical(c_bool), parameter :: true=.true., false=.false.
    integer(local_ordinal_type) :: lclrow, cols(1), numvalid
    integer(global_ordinal_type) :: gblrow, gblcols(1)
    real(scalar_type), parameter :: zero=0.
    real(scalar_type) :: vals(1), zeros(1)=[zero]

    OUT0("Starting TpetraCrsMatrix_ActiveFill")

    ! TODO: Creating the matrix with the two map constructor is causing some problems.  Could just be my implementation
    ! create Map
    call Map%create(invalid, ione, index_base, comm); TEST_IERR()

    call Mat%create(map, map, izero, DynamicProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    lclrow = 1; cols(1) = 1; vals(1) = 0.
    !TODO: In all of the following, use *Local* methods, not Global
    ! call Mat%insertLocalValues(lclrow, cols, tuple<Scalar>(0)); TEST_IERR() TODO: Reinstate
    gblrow = map%getGlobalElement(lclrow)  ! TODO: Remove
    gblcols(1) = gblrow  ! TODO: Remove
    call Mat%insertGlobalValues(gblrow, gblcols, zeros); TEST_IERR()  ! TODO: Remove

    call params%create("ANONOMOUS")
    ! call params%set("Optimize Storage", false) TODO: boolean parameters
    call Mat%fillComplete(params);
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    numvalid = 0

    ! It's forbidden to call any of the *LocalValues methods if the
    ! matrix is fill complete (not fill active).

    ! TEST_THROW(call Mat%insertLocalValues(lclrow, cols, vals)) TODO: Reinstate
    ! TODO: The following line should set ierr/=0, but it just throws from Tpetra
    ! TEST_THROW(call Mat%insertGlobalValues(gblrow, gblcols, vals)) ! TODO: Remove

    ! numvalid = Mat%replaceLocalValues(lclrow, cols, vals); TEST_IERR() ! TODO: Reinstate
    numvalid = Mat%replaceGlobalValues(gblrow, gblcols, vals); TEST_IERR() ! TODO: Remove
    TEST_EQUALITY(numvalid, invalid)

    ! numvalid = Mat%sumIntoLocalValues(lcrow, cols, vals); TEST_IERR() ! TODO: Reinstate
    numvalid = Mat%sumIntoGlobalValues(gblrow, gblcols, vals); TEST_IERR() ! TODO: Remove
    TEST_EQUALITY(numvalid, invalid)

    ! TODO: The following should set ierr/=0, but it just throws from Tpetra
    !TEST_THROW(call Mat%setAllToScalar(zero))
    !TEST_THROW(call Mat%scale(zero))
    !TEST_THROW(call Mat%globalAssemble())
    !TEST_THROW(call Mat%fillComplete())

    call params%release()
    call Mat%release()

    call Mat%create(map, map, izero, DynamicProfile); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    lclrow = 1; cols(1) = 1; vals(1) = zero;
    !call Mat%insertLocalValues(lclrow, cols, vals); TEST_IERR() TODO: Reinstate
    call Mat%insertGlobalValues(gblrow, gblcols, vals); TEST_IERR() ! TODO: Remove

    call params%create("ANONOMOUS"); TEST_IERR()
    !call params%set("Optimize Storage", .false.); TEST_IERR() TODO: boolean parameters
    call Mat%fillComplete(params); TEST_IERR()
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call Mat%resumeFill(); TEST_IERR()
    TEST_ASSERT(Mat%isFillActive())
    TEST_ASSERT((.not. Mat%isFillComplete()))
    !TEST_NOTHROW(call Mat%insertLocalValues(lclrow, cols, vals)) TODO: Reinstate
    !TEST_NOTHROW(numvalid = Mat%replaceLocalValues(lclrow, cols, vals)) TODO: Reinstate
    !TEST_NOTHROW(call Mat%sumIntoLocalValues(lclrow, cols, vals)) TODO: Reinstate

    ! TODO: The following should NOT set ierr/=0 but does
    !TEST_NOTHROW(call Mat%insertGlobalValues(gblrow, gblcols, vals)) ! TODO: Remove
    TEST_NOTHROW(numvalid = Mat%replaceGlobalValues(gblrow, gblcols, vals)) ! TODO: Remove
    TEST_NOTHROW(numvalid = Mat%sumIntoGlobalValues(gblrow, gblcols, vals)) ! TODO: Remove
    TEST_NOTHROW(call Mat%setAllToScalar(zero))
    TEST_NOTHROW(call Mat%scale(zero))
    TEST_NOTHROW(call Mat%globalAssemble())

    TEST_NOTHROW(call Mat%fillComplete())
    TEST_ASSERT((.not. Mat%isFillActive()))
    TEST_ASSERT(Mat%isFillComplete())

    call params%release()
    call Map%release()
    call Mat%release()

    OUT0("Finished TpetraCrsMatrix_ActiveFill!")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_ActiveFill)

  ! -------------------------------getCrsGraph-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getCrsGraph)
    type(TpetraCrsMatrix) :: Obj
    OUT0("Starting TpetraCrsMatrix_getCrsGraph")

    success = .false.

    !call Obj%create(); TEST_IERR()
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

    !call Obj%create(); TEST_IERR()
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

    !call Obj%create(); TEST_IERR()
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

    !call Obj%create(); TEST_IERR()
    !fresult = Obj%supportsRowViews(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_supportsRowViews: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_supportsRowViews")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_supportsRowViews)

  ! -----------------------------getGlobalRowView----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalRowView)
    type(TpetraCrsMatrix) :: Obj
    integer(C_LONG_LONG) :: globalrow
    type(TeuchosArrayViewLongLongConst) :: indices
    type(TeuchosArrayViewDoubleConst) :: values
    OUT0("Starting TpetraCrsMatrix_getGlobalRowView")

    success = .false.

    globalrow = 0
    !call indices%create(); TEST_IERR()
    !call values%create(); TEST_IERR()
    !call Obj%create(); TEST_IERR()
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
    integer(C_INT) :: lclrow
    integer(C_INT) :: nument
    type(C_PTR) :: lclcolinds
    type(C_PTR) :: vals
    OUT0("Starting TpetraCrsMatrix_getLocalRowViewRaw")

    success = .false.

    lclrow = 0
    nument = 0
    !call lclcolinds%create(); TEST_IERR()
    !call vals%create(); TEST_IERR()
    !call Obj%create(); TEST_IERR()
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
    real(C_DOUBLE) :: dampingfactor
    integer(kind(ESweepDirection)) :: direction
    integer(C_INT) :: numsweeps
    OUT0("Starting TpetraCrsMatrix_gaussSeidel")

    success = .false.

    !call b%create(); TEST_IERR()
    !call x%create(); TEST_IERR()
    !call d%create(); TEST_IERR()
    dampingfactor = 0
    direction = 0
    numsweeps = 0
    !call Obj%create(); TEST_IERR()
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
    type(TeuchosArrayViewInt) :: rowindices
    real(C_DOUBLE) :: dampingfactor
    integer(kind(ESweepDirection)) :: direction
    integer(C_INT) :: numsweeps
    OUT0("Starting TpetraCrsMatrix_reorderedGaussSeidel")

    success = .false.

    !call b%create(); TEST_IERR()
    !call x%create(); TEST_IERR()
    !call d%create(); TEST_IERR()
    !call rowindices%create(); TEST_IERR()
    dampingfactor = 0
    direction = 0
    numsweeps = 0
    !call Obj%create(); TEST_IERR()
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
    real(C_DOUBLE) :: dampingfactor
    integer(kind(ESweepDirection)) :: direction
    integer(C_INT) :: numsweeps
    logical(C_BOOL) :: zeroinitialguess
    OUT0("Starting TpetraCrsMatrix_gaussSeidelCopy")

    success = .false.

    !call x%create(); TEST_IERR()
    !call b%create(); TEST_IERR()
    !call d%create(); TEST_IERR()
    dampingfactor = 0
    direction = 0
    numsweeps = 0
    zeroinitialguess = .false.
    !call Obj%create(); TEST_IERR()
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
    type(TeuchosArrayViewInt) :: rowindices
    real(C_DOUBLE) :: dampingfactor
    integer(kind(ESweepDirection)) :: direction
    integer(C_INT) :: numsweeps
    logical(C_BOOL) :: zeroinitialguess
    OUT0("Starting TpetraCrsMatrix_reorderedGaussSeidelCopy")

    success = .false.

    !call x%create(); TEST_IERR()
    !call b%create(); TEST_IERR()
    !call d%create(); TEST_IERR()
    !call rowindices%create(); TEST_IERR()
    dampingfactor = 0
    direction = 0
    numsweeps = 0
    zeroinitialguess = .false.
    !call Obj%create(); TEST_IERR()
    !call Obj%reorderedGaussSeidelCopy(x, b, d, rowindices, dampingfactor, direction, numsweeps, zeroinitialguess); TEST_IERR()

    !call x%release(); TEST_IERR()
    !call b%release(); TEST_IERR()
    !call d%release(); TEST_IERR()
    !call rowindices%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_reorderedGaussSeidelCopy: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_reorderedGaussSeidelCopy")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_reorderedGaussSeidelCopy)

  ! -----------------------------getGlobalRowCopy----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalRowCopy)
    type(TpetraCrsMatrix) :: Obj
    integer(C_LONG_LONG) :: globalrow
    integer(C_LONG_LONG), allocatable :: indices(:)
    real(C_DOUBLE), allocatable :: values(:)
    integer(C_SIZE_T) :: numindices
    OUT0("Starting TpetraCrsMatrix_getGlobalRowCopy")

    success = .false.

    globalrow = 0
    !allocate(indices(:)(0))
    !allocate(values(:)(0))
    numindices = 0
    !call Obj%create(); TEST_IERR()
    !call Obj%getGlobalRowCopy(globalrow, indices(:), values(:), numindices); TEST_IERR()

    !deallocate(indices(:))
    !deallocate(values(:))
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_getGlobalRowCopy: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_getGlobalRowCopy")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_getGlobalRowCopy)

  ! ------------------------------replaceColMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_replaceColMap)
    type(TpetraCrsMatrix) :: Obj
    type(TpetraMap) :: newcolmap
    OUT0("Starting TpetraCrsMatrix_replaceColMap")

    success = .false.

    !call newcolmap%create(); TEST_IERR()
    !call Obj%create(); TEST_IERR()
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

    !call newdomainmap%create(); TEST_IERR()
    !call newimporter%create(); TEST_IERR()
    !call Obj%create(); TEST_IERR()
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

    !call newmap%create(); TEST_IERR()
    !call Obj%create(); TEST_IERR()
    !call Obj%removeEmptyProcessesInPlace(newmap); TEST_IERR()

    !call newmap%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraCrsMatrix_removeEmptyProcessesInPlace: Test not yet implemented'

    OUT0("Finished TpetraCrsMatrix_removeEmptyProcessesInPlace")

  END_FORTRILINOS_UNIT_TEST(TpetraCrsMatrix_removeEmptyProcessesInPlace)


end program test_TpetraCrsMatrix
