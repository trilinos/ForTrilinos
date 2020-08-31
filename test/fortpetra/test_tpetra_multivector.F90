! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program test_TpetraMultiVector
#include "ForTrilinos_config.h"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra
  use test_Tpetra_multivector_helper

  implicit none
  type(TeuchosComm) :: comm
  integer, parameter :: dp = kind(0.d0)
  character(len=256), parameter :: FILENAME="test_tpetra_multivector.f90"

  SETUP_TEST()

#if FORTRILINOS_USE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm(); FORTRILINOS_CHECK_IERR()
#endif

  ADD_SUBTEST_AND_RUN(TpetraMultiVector_basic)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_zeroScaleUpdate)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_countNormInf)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_norm2)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_replaceMap)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_reciprocal)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_abs)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_description)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_meanValue)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_multiply)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_reduce)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_replaceGlobalValue)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_replaceLocalValue)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_get1dCopy)

  ADD_SUBTEST_AND_RUN(TpetraMultiVector_swap)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_putScalar)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_subCopy)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_subView)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_subViewNonConst)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_offsetView)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_offsetViewNonConst)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_getData)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_getDataNonConst)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_get1dView)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_get1dViewNonConst)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_getNumVectors)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_getLocalLength)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_getGlobalLength)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_getStride)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_isConstantStride)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_isSameSize)

  ! These tests are Cuda so will develop them on the new
  ! multi-node branch and merge it all later.
  ! ADD_SUBTEST_AND_RUN(TpetraMultiVector_sync_host)
  ! ADD_SUBTEST_AND_RUN(TpetraMultiVector_sync_device)
  ! ADD_SUBTEST_AND_RUN(TpetraMultiVector_need_sync_host)
  ! ADD_SUBTEST_AND_RUN(TpetraMultiVector_need_sync_device)
  ! ADD_SUBTEST_AND_RUN(TpetraMultiVector_modify_device)
  ! ADD_SUBTEST_AND_RUN(TpetraMultiVector_modify_host)


  ! The following methods are deprecated
  !ADD_SUBTEST_AND_RUN(TpetraMultiVector_normWeighted)

  ! The following methods are listed as "may change or disappear at any time"
  !ADD_SUBTEST_AND_RUN(TpetraMultiVector_setCopyOrView)
  !ADD_SUBTEST_AND_RUN(TpetraMultiVector_getCopyOrView)
  !ADD_SUBTEST_AND_RUN(TpetraMultiVector_removeEmptyProcessesInPlace)

  call comm%release();  TEST_IERR()

  TEARDOWN_TEST()

contains

  ! --------------------------------Basic------------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Basic)
    Type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    integer(size_type), parameter :: num_vecs=12
    integer, parameter :: num_local=2
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    Vec = TpetraMultiVector(map, num_vecs); TEST_IERR()

    TEST_EQUALITY(Vec%getNumVectors(), num_vecs)
    TEST_EQUALITY(Vec%getLocalLength(), num_local)
    TEST_EQUALITY(Vec%getGlobalLength(),
                  int(num_local*comm%getSize(), size_type))
    TEST_EQUALITY(Vec%getStride(), int(num_local, size_type))
    TEST_ASSERT(Vec%isConstantStride())

    call Vec%release(); TEST_IERR()
    call map%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Basic)

  ! -----------------------------ZeroScaleUpdate------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_ZeroScaleUpdate)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: A, B, A2, C
    integer(size_type), parameter :: num_vecs=2, LDA=2
    integer, parameter :: num_local=2
    logical :: zeroout
    real(scalar_type) :: norms(num_vecs), zeros(num_vecs), values(6)

    zeros = 0.d0

    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()

    ! values = {1, 1, 2, 2, 4, 4}
    ! values(1:4) = {1, 1, 2, 2} = [1 2]
    !                            = [1 2]
    ! values(3:)  = {2, 2, 4, 4} = [2 4]
    !                            = [2 4]
    ! a multivector A constructed from the first
    ! has values .5 of a multivector B constructed from the second
    ! then 2*A - B = 0
    ! we test both scale(), both update(), and norm()
    values = [1.d0, 1.d0, 2.d0, 2.d0, 4.d0, 4.d0]

    ! TODO: Multivec create to take array, not ArrayView
    A = TpetraMultiVector(map, values(1:4), LDA, num_vecs); TEST_IERR()
    B = TpetraMultiVector(map, values(3:), LDA, num_vecs); TEST_IERR()

    !
    !      [.... ....]
    ! A == [ones ones]
    !      [.... ....]
    !
    !      [.... ....]
    ! B == [twos twos]
    !      [.... ....]
    !
    !   set A2 = A
    !   scale it by 2 in situ
    !   check that it equals B: subtraction in situ
    A2 = TpetraMultiVector(A, TeuchosCopy); TEST_IERR()
    call A2%scale(2.d0)
    call A2%update(-1.d0, B, 1.0d0)
    call A2%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, zeros, epsilon(0.d0))
    call A2%release();  TEST_IERR()

    ! set A2 = A
    ! check that it equals B: scale, subtraction in situ
    A2 = TpetraMultiVector(A, TeuchosCopy); TEST_IERR()
    call A2%update(-1.d0, B, 2.d0)
    call A2%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, zeros, epsilon(0.d0))
    call A2%release();  TEST_IERR()

    ! set C random
    ! set it to zero by combination with A,B
    zeroout = .false.
    C = TpetraMultiVector(map, num_vecs, zeroout); TEST_IERR()
    call C%randomize()
    call C%update(-1.d0, B, 2.d0, A, 0.d0)
    call C%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, zeros, epsilon(0.d0))
    call C%release();  TEST_IERR()

    ! set C random
    ! scale it ex-situ
    ! check that it equals B: subtraction in situ
    C = TpetraMultiVector(map, num_vecs, zeroout); TEST_IERR()
    call C%scale(2.d0, A)
    call C%update(1.0d0, B, -1.d0)
    call C%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, zeros, epsilon(0.d0))
    call C%release();  TEST_IERR()

    ! Clean up
    call A%release(); TEST_IERR()
    call B%release(); TEST_IERR()
    call map%release(); TEST_IERR()

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_ZeroScaleUpdate)

  ! ----------------------------CountNormInf---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_CountNormInf)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    integer(size_type), parameter :: num_vecs=3, LDA=2
    integer, parameter :: num_local=2
    real(scalar_type) :: values(num_vecs*num_local), answer(num_vecs), norms(num_vecs)

    OUT0("Starting CountNormInf")

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    ! values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    !                               [0 1 2]
    ! normInf(values) = [0 1 2]
    ! over all procs, this is [0 1 2]
    values(:) = [0., 0., 1., 1., 2., 2.]
    Vec = TpetraMultiVector(map, values, LDA, num_vecs); TEST_IERR()
    answer(:) = [0., 1., 2.]

    ! do the dots
    call Vec%normInf(norms)

    call Vec%release();  TEST_IERR()
    call map%release();  TEST_IERR()

    ! check the answers
    TEST_FLOATING_ARRAY_EQUALITY(norms, answer, epsilon(answer(1)))

    OUT0("Finished CountNormInf!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_CountNormInf)

  ! --------------------------------Norm2------------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Norm2)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    integer(size_type), parameter :: num_vecs=7
    integer, parameter :: num_local=13
    real(scalar_type) :: norms_rand(num_vecs), norms_zero(num_vecs)

    OUT0("Starting Norm2")

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()

    Vec = TpetraMultiVector(map, num_vecs); TEST_IERR()
    call Vec%randomize(); TEST_IERR()

    ! Take the norms, they should not be zero
    call Vec%norm2(norms_rand)

    ! Zero the vector
    call Vec%putScalar(0.d0)

    ! Take the norms, they should be zero
    call Vec%norm2(norms_zero)

    ! Check the answers
    TEST_FLOATING_ARRAY_INEQUALITY(norms_rand, 0.d0, epsilon(0.d0))
    TEST_FLOATING_ARRAY_EQUALITY(norms_zero, 0.d0, epsilon(0.d0))

    call Vec%release();  TEST_IERR()
    call map%release();  TEST_IERR()

    OUT0("Finished Norm2!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Norm2)

  ! -----------------------------Reciprocal----------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Reciprocal)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: A, B
    integer(size_type), parameter :: num_vecs=2
    integer, parameter :: num_local=10
    real(scalar_type) :: dots(num_vecs)

    OUT0("Starting Reciprocal")

    ! create a Map
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()

    A = TpetraMultiVector(map, num_vecs); TEST_IERR()
    B = TpetraMultiVector(map, num_vecs); TEST_IERR()
    call A%putScalar(5.0d0); TEST_IERR()
    call B%reciprocal(A); TEST_IERR()

    ! Take the dots, they should 1.0d0
    call A%dot(B, dots)

    ! Check the answers
    TEST_FLOATING_ARRAY_EQUALITY(dots, 1.0d0*num_local*comm%getSize(), epsilon(1.0d0))

    call A%release();  TEST_IERR()
    call B%release();  TEST_IERR()
    call map%release();  TEST_IERR()

    OUT0("Finished Reciprocal!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Reciprocal)

  ! --------------------------------ReplaceMap-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_ReplaceMap)
    type(TpetraMultiVector) :: Obj
    type(TpetraMap) :: map1, map2
    integer(size_type), parameter :: num_vecs=1
    integer(global_ordinal_type) :: num_global

    OUT0("Starting ReplaceMap")
    ! NOTE: This test only tests the interface does not throw.
    ! It does not test correctness
    num_global = 4 * comm%getSize()
    map1 = TpetraMap(num_global, comm); TEST_IERR()

    num_global = 5 * comm%getSize()
    map2 = TpetraMap(num_global, comm); TEST_IERR()

    Obj = TpetraMultiVector(map1, num_vecs); TEST_IERR()

    call Obj%replaceMap(map2); TEST_IERR()

    call map1%release(); TEST_IERR()
    call map2%release(); TEST_IERR()

    call Obj%release(); TEST_IERR()

    OUT0("Finished ReplaceMap")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_ReplaceMap)

  ! ----------------------------------Abs------------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Abs)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: A, B, A2
    integer(size_type), parameter :: num_vecs=2
    integer, parameter :: num_local=10
    real(scalar_type) :: norms(num_vecs)

    OUT0("Starting Abs")
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()

    A = TpetraMultiVector(map, num_vecs); TEST_IERR()
    B = TpetraMultiVector(map, num_vecs); TEST_IERR()
    call A%putScalar(-1.d0); TEST_IERR()
    call A%abs(B)

    !   set A2 = A
    !   scale it by 2 in situ
    !   check that it equals B: subtraction in situ
    A2 = TpetraMultiVector(A, TeuchosCopy); TEST_IERR()
    call A2%update(1.0d0, B, 1.0d0)
    call A2%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, 0.d0, epsilon(0.d0))

    call A2%release();  TEST_IERR()
    call B%release();  TEST_IERR()
    call A%release();  TEST_IERR()
    call map%release();  TEST_IERR()

    OUT0("Finished Abs")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Abs)

! --------------------------------description--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Description)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    character(kind=C_CHAR, len=:), allocatable :: description
    integer(size_type), parameter :: num_vecs=2
    integer, parameter :: num_local=10
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    Vec = TpetraMultiVector(map, num_vecs)
    description = Vec%description(); TEST_IERR()
    write(*,*) "Vector description:", description
    call Vec%release(); TEST_IERR()
    call map%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Description)

  ! --------------------------------MeanValue--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_MeanValue)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    integer(size_type), parameter :: num_vecs=2, LDA=2
    integer, parameter :: num_local=2
    real(scalar_type) :: values(4), means(num_vecs), answer(num_vecs)

    OUT0("Starting MeanValue")

    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()

    ! values = {2, 6, 3, 1} = [2 3]
    !                         [6 1]
    values(:) = [2., 6., 3., 1.]

    Vec = TpetraMultiVector(map, values, LDA, num_vecs); TEST_IERR()

    !means = xxx(); TEST_IERR()
    call Vec%meanValue(means); TEST_IERR()

    answer = [4., 2.]
    TEST_FLOATING_ARRAY_EQUALITY(means, answer, epsilon(means(1)))

    call Vec%release(); TEST_IERR()
    call map%release(); TEST_IERR()

    OUT0("Finished MeanValue")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_MeanValue)

  ! ---------------------------------Multiply--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Multiply)
    type(TpetraMap) :: map2n, map3n, lmap2, lmap3
    type(TpetraMultiVector) :: mv3nx2, mv3nx3, mv2x2, mv2x3, mv3x2, mv3x3
    real(scalar_type), parameter :: S0=0., S1=1.
    integer(int_type) :: num_images
    real(scalar_type) :: check(9)

    OUT0("Starting Multiply")

    map2n = TpetraMap(TPETRA_GLOBAL_INVALID, 2, comm); TEST_IERR()
    map3n = TpetraMap(TPETRA_GLOBAL_INVALID, 3, comm); TEST_IERR()

    lmap2 = TpetraMap(2_size_type, comm, TpetraLocallyReplicated); TEST_IERR()
    lmap3 = TpetraMap(3_size_type, comm, TpetraLocallyReplicated); TEST_IERR()

    mv3nx2 = TpetraMultiVector(map3n, 2_size_type)
    mv3nx3 = TpetraMultiVector(map3n, 3_size_type)

    mv2x2 = TpetraMultiVector(lmap2, 2_size_type)
    mv2x3 = TpetraMultiVector(lmap2, 3_size_type)
    mv3x2 = TpetraMultiVector(lmap3, 2_size_type)
    mv3x3 = TpetraMultiVector(lmap3, 3_size_type)

    num_images = comm%getSize()
    check = 3 * num_images

    call mv2x2%multiply(TEUCHOSCONJ_TRANS,TEUCHOSNO_TRANS,S1,mv3nx2,mv3nx2,S0); TEST_IERR()
    call mv2x3%multiply(TEUCHOSCONJ_TRANS,TEUCHOSNO_TRANS,S1,mv3nx2,mv3nx3,S0); TEST_IERR()
    call mv3x2%multiply(TEUCHOSCONJ_TRANS,TEUCHOSNO_TRANS,S1,mv3nx3,mv3nx2,S0); TEST_IERR()
    call mv3x3%multiply(TEUCHOSCONJ_TRANS,TEUCHOSNO_TRANS,S1,mv3nx3,mv3nx3,S0); TEST_IERR()

    !a = fortran array view of the data
    !n = size of a
    !TEST_FLOATING_ARRAY_EQUALITY(a, check(1:n), epsilon(a(1)))

    call map2n%release(); TEST_IERR()
    call map3n%release(); TEST_IERR()
    call lmap2%release(); TEST_IERR()
    call lmap3%release(); TEST_IERR()
    call mv3nx2%release(); TEST_IERR()
    call mv3nx3%release(); TEST_IERR()
    call mv2x2%release(); TEST_IERR()
    call mv2x3%release(); TEST_IERR()
    call mv3x2%release(); TEST_IERR()
    call mv3x3%release(); TEST_IERR()

    OUT0("Finished Multiply!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Multiply)

  ! ----------------------------------Reduce---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Reduce)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    integer(global_size_type) :: num_global
    integer(size_type), parameter :: num_vecs=1
    integer, parameter :: num_local=2
    OUT0("Starting Reduce")
    num_global = num_local * comm%getSize()
    map = TpetraMap(num_global, comm, TpetraLocallyReplicated); TEST_IERR()
    Vec = TpetraMultiVector(map, num_vecs); TEST_IERR()
    call Vec%putScalar(2.d0); TEST_IERR()
    call Vec%reduce(); TEST_IERR()
    call Vec%release(); TEST_IERR()
    call map%release(); TEST_IERR()
    OUT0("Finished Reduce!")
  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Reduce)

  ! ----------------------------replaceGlobalValue---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_ReplaceGlobalValue)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec, OneV
    integer(size_type) :: col
    integer(size_type), parameter :: num_vecs=2
    integer, parameter :: num_local=4
    integer :: lclrow
    real(scalar_type) :: value, expected, dots(num_vecs)
    integer(global_ordinal_type) :: gblrow, num_global
    OUT0("Starting replaceGlobalValue")
    num_global = num_local * comm%getSize()
    map = TpetraMap(num_global, comm); TEST_IERR()

    Vec = TpetraMultiVector(map, num_vecs); TEST_IERR()
    OneV = TpetraMultiVector(map, num_vecs); TEST_IERR()
    call OneV%putScalar(1.0d0); TEST_IERR()

    do lclrow = 1, num_local
      gblrow = map%getGlobalElement(lclrow)
      value = real(gblrow, kind=scalar_type)
      do col = 1, num_vecs
        call Vec%replaceGlobalValue(gblrow, col, value); TEST_IERR()
      end do
    end do

    call Vec%dot(OneV, dots)
    expected = real(num_global * (num_global + 1) / 2, kind=scalar_type)
    TEST_FLOATING_EQUALITY(dots(1), dots(2), epsilon(dots(2)))
    TEST_FLOATING_EQUALITY(expected, dots(1), epsilon(dots(2)))

    call Vec%release(); TEST_IERR()
    call OneV%release(); TEST_IERR()
    call map%release(); TEST_IERR()

    OUT0("Finished replaceGlobalValue!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_ReplaceGlobalValue)

  ! ----------------------------ReplaceLocalValue----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_ReplaceLocalValue)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec, OneV
    integer(size_type) :: col
    integer :: lclrow
    integer(size_type), parameter :: num_vecs=2
    integer, parameter :: num_local=4
    real(scalar_type) :: value, dots(num_vecs), expected
    OUT0("Starting ReplaceLocalValue")

    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    Vec = TpetraMultiVector(map, num_vecs); TEST_IERR()
    OneV = TpetraMultiVector(map, num_vecs); TEST_IERR()
    call OneV%putScalar(1.0d0)

    do lclrow = 1, num_local
      value = real(lclrow, kind=scalar_type)
      do col = 1, num_vecs
        call Vec%replaceLocalValue(lclrow, col, value); TEST_IERR()
      end do
    end do

    call Vec%dot(OneV, dots)
    expected = real(comm%getSize() * (num_local * (num_local + 1)) / 2, kind=scalar_type)
    TEST_FLOATING_EQUALITY(dots(1), dots(2), epsilon(dots(2)))
    TEST_FLOATING_EQUALITY(expected, dots(1), epsilon(dots(2)))

    call Vec%release(); TEST_IERR()
    call OneV%release(); TEST_IERR()
    call map%release();  TEST_IERR()

    OUT0("Finished ReplaceLocalValue!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_ReplaceLocalValue)

  ! --------------------------------Get1dCopy--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Get1dCopy)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    integer(size_type) :: lda
    real(scalar_type), allocatable :: a(:)
    integer(size_type), parameter :: num_vecs=2
    integer, parameter :: num_local=4
    OUT0("Starting Get1dCopy")
    map = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm); TEST_IERR()
    Vec = TpetraMultiVector(map, num_vecs); TEST_IERR()
    call Vec%putScalar(1.0d0); TEST_IERR()
    allocate(a(num_vecs*num_local*comm%getSize()))
    a = 0.
    lda = num_local*comm%getSize()
    call Vec%Get1dCopy(a, lda); TEST_IERR()
    call Vec%release(); TEST_IERR()

    TEST_FLOATING_ARRAY_EQUALITY(a, 1.0d0, epsilon(1.0d0))
    deallocate(a)

    call map%release(); TEST_IERR()
    call Vec%release(); TEST_IERR()

    OUT0("Finished Get1dCopy!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Get1dCopy)

  ! -----------------------------------swap----------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_swap)
    type(TpetraMultiVector) :: A, B
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2
    real(scalar_type) :: normsA(num_vecs), normsB(num_vecs)

    OUT0("Starting TpetraMultiVector_swap!")

    ! Create test vectors, fill one of them with 1.0
    call Tpetra_MV_Create(comm, num_local, num_vecs, A)
    call Tpetra_MV_Create(comm, num_local, num_vecs, B, 1.0_dp)

    ! Check the entries are what they should be
    call A%norm2(normsA)
    call B%norm2(normsB)
    TEST_FLOATING_ARRAY_EQUALITY(normsA, 0.0_dp, epsilon(0.0_dp))
    TEST_FLOATING_ARRAY_INEQUALITY(normsB, 0.0_dp, epsilon(0.0_dp))

    ! Swap entries and ensure it went through
    call A%swap(B)
    call A%norm2(normsA)
    call B%norm2(normsB)
    TEST_FLOATING_ARRAY_INEQUALITY(normsA, 0.0_dp, epsilon(0.0_dp))
    TEST_FLOATING_ARRAY_EQUALITY(normsB, 0.0_dp, epsilon(0.0_dp))

    call A%release(); TEST_IERR()
    call B%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_swap!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_swap)

  ! --------------------------------putScalar--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_putScalar)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2
    real(dp) :: val, expect
    real(scalar_type) :: mynorms(num_vecs)

    OUT0("Starting TpetraMultiVector_putScalar!")

    val = 1.0_dp
    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    call vec%norm1(mynorms)
    TEST_FLOATING_ARRAY_EQUALITY(mynorms, 0.0_dp, epsilon(0.0_dp))

    call vec%putScalar(val)
    call vec%norm1(mynorms)
    expect = 1.0_dp * real(num_local, kind = dp) * comm%getSize()
    TEST_FLOATING_ARRAY_EQUALITY(mynorms, expect, epsilon(expect))

    call vec%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_putScalar!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_putScalar)

  ! ---------------------------------subCopy---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_subCopy)
    type(TpetraMultiVector) :: vec, res
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2
    integer(size_type) :: cols(1)
    real(dp) :: resNorm(1), expect

    OUT0("Starting TpetraMultiVector_subCopy!")

    cols(1) = 1
    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    res = vec%subCopy(cols)

    TEST_EQUALITY(res%getNumVectors(), 1_size_type)
    TEST_EQUALITY(res%getLocalLength(), num_local)
    call res%norm1(resNorm)
    TEST_FLOATING_ARRAY_EQUALITY(resNorm, 0.0_dp, epsilon(0.0_dp))

    call vec%putScalar(1.0_dp)
    call res%norm1(resNorm)
    ! vec's data changed, but res is its own vector, so no change to norm
    TEST_FLOATING_ARRAY_EQUALITY(resNorm, 0.0_dp, epsilon(0.0_dp))

    call vec%release(); TEST_IERR()
    call res%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_subCopy!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_subCopy)

  ! ---------------------------------subView---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_subView)
    type(TpetraMultiVector) :: vec, res
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2
    integer(size_type) :: cols(1)
    real(dp) :: resNorm(1), vecNorm(num_vecs), expect

    OUT0("Starting TpetraMultiVector_subView!")

    cols(1) = 1
    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    res = vec%subView(cols)

    TEST_EQUALITY(res%getNumVectors(), 1_size_type)
    TEST_EQUALITY(res%getLocalLength(), num_local)
    call res%norm1(resNorm)
    TEST_FLOATING_ARRAY_EQUALITY(resNorm, 0.0_dp, epsilon(0.0_dp))

    call vec%putScalar(1.0_dp)
    expect = 1.0_dp * num_local * comm%getSize()
    call res%norm1(resNorm)
    ! vec's data changed, so res should change with it
    TEST_FLOATING_ARRAY_EQUALITY(resNorm, expect, epsilon(expect))

    call res%putScalar(0.0_dp)
    call vec%norm1(vecNorm)
    TEST_FLOATING_EQUALITY(vecNorm(1), 0.0_dp, epsilon(0.0_dp))
    TEST_FLOATING_EQUALITY(vecNorm(2), expect, epsilon(0.0_dp))

    call vec%release(); TEST_IERR()
    call res%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_subView!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_subView)

  ! -----------------------------subViewNonConst------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_subViewNonConst)
    type(TpetraMultiVector) :: vec, res
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2
    integer(size_type) :: cols(1)
    real(dp) :: resNorm(1), vecNorm(num_vecs), expect

    OUT0("Starting TpetraMultiVector_subViewNonConst!")

    cols(1) = 1
    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    res = vec%subViewNonConst(cols)

    TEST_EQUALITY(res%getNumVectors(), 1_size_type)
    TEST_EQUALITY(res%getLocalLength(), num_local)
    call res%norm1(resNorm)

    TEST_FLOATING_ARRAY_EQUALITY(resNorm, 0.0_dp, epsilon(0.0_dp))

    call vec%putScalar(1.0_dp)
    expect = 1.0_dp * num_local * comm%getSize()
    call res%norm1(resNorm)
    ! vec's data changed, so res should change with it
    TEST_FLOATING_ARRAY_EQUALITY(resNorm, expect, epsilon(expect))

    call res%putScalar(0.0_dp)
    call vec%norm1(vecNorm)
    TEST_FLOATING_EQUALITY(vecNorm(1), 0.0_dp, epsilon(0.0_dp))
    TEST_FLOATING_EQUALITY(vecNorm(2), expect, epsilon(0.0_dp))

    call vec%release(); TEST_IERR()
    call res%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_subViewNonConst!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_subViewNonConst)

  ! --------------------------------offsetView-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_offsetView)
    type(TpetraMultiVector) :: vec, res1, res2
    type(TpetraMap) :: split
    integer, parameter :: num_local = 4, num_split = 2
    integer(size_type), parameter :: num_vecs = 2
    integer(size_type) :: offset
    real(dp) :: res1Norm(2), res2Norm(2), expect

    OUT0("Starting TpetraMultiVector_offsetView!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    split = TpetraMap(TPETRA_GLOBAL_INVALID, num_split, comm)

    res1 = vec%offsetView(split, 0_size_type)
    res2 = vec%offsetView(split, int(res1%getLocalLength(), kind = size_type))

    TEST_EQUALITY(res1%getLocalLength() + res2%getLocalLength(), vec%getLocalLength())
    TEST_EQUALITY(res1%getNumVectors(), vec%getNumVectors())
    TEST_EQUALITY(res2%getNumVectors(), vec%getNumVectors())

    call vec%putScalar(1.0_dp)
    call res1%norm1(res1Norm)
    call res2%norm1(res2Norm)

    expect = 1.0_dp * num_split * comm%getSize()
    TEST_FLOATING_ARRAY_EQUALITY(res1Norm, expect, epsilon(expect))
    TEST_FLOATING_ARRAY_EQUALITY(res2Norm, expect, epsilon(expect))

    call vec%release(); TEST_IERR()
    call res1%release(); TEST_IERR()
    call res2%release(); TEST_IERR()
    call split%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_offsetView!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_offsetView)

  ! ----------------------------offsetViewNonConst---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_offsetViewNonConst)

    type(TpetraMultiVector) :: vec, res1, res2
    type(TpetraMap) :: split
    integer, parameter :: num_local = 4, num_split = 2
    integer(size_type), parameter :: num_vecs = 2
    integer(size_type) :: offset
    real(dp) :: res1Norm(2), res2Norm(2), expect

    OUT0("Starting TpetraMultiVector_offsetViewNonConst!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    split = TpetraMap(TPETRA_GLOBAL_INVALID, num_split, comm)

    res1 = vec%offsetView(split, 0_size_type)
    res2 = vec%offsetView(split, int(res1%getLocalLength(), kind = size_type))

    TEST_EQUALITY(res1%getLocalLength() + res2%getLocalLength(), vec%getLocalLength())
    TEST_EQUALITY(res1%getNumVectors(), vec%getNumVectors())
    TEST_EQUALITY(res2%getNumVectors(), vec%getNumVectors())

    call vec%putScalar(1.0_dp)
    call res1%norm1(res1Norm)
    call res2%norm1(res2Norm)

    expect = 1.0_dp * num_split * comm%getSize()
    TEST_FLOATING_ARRAY_EQUALITY(res1Norm, expect, epsilon(expect))
    TEST_FLOATING_ARRAY_EQUALITY(res2Norm, expect, epsilon(expect))

    call vec%release(); TEST_IERR()
    call res1%release(); TEST_IERR()
    call res2%release(); TEST_IERR()
    call split%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_offsetViewNonConst!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_offsetViewNonConst)

  ! ---------------------------------getData---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_getData)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2
    real(dp) :: a(num_local), vecNorm(num_vecs), expect(num_vecs)
    integer(size_type) :: i, j

    OUT0("Starting TpetraMultiVector_getData!")

    j = 1
    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    a = vec%getData(j)

    a = 1.0_dp

    call vec%norm1(vecNorm)
    expect = 0.0_dp
    expect(1) = 1.0_dp * num_local * comm%getSize()

    TEST_FLOATING_ARRAY_EQUALITY(vecNorm, expect, epsilon(expect(1)))

    call vec%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_getData!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_getData)

  ! -----------------------------getDataNonConst------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_getDataNonConst)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2
    real(dp) :: a(num_local), vecNorm(num_vecs), expect(num_vecs)
    integer(size_type) :: j

    OUT0("Starting TpetraMultiVector_getDataNonConst!")

    j = 1
    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    a = vec%getDataNonConst(j)

    a = 1.0_dp

    call vec%norm1(vecNorm)
    expect = 0.0_dp
    expect(1) = 1.0_dp * num_local * comm%getSize()

    TEST_FLOATING_ARRAY_EQUALITY(vecNorm, expect, epsilon(expect(1)))

    call vec%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_getDataNonConst!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_getDataNonConst)

  ! --------------------------------get1dView--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_get1dView)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2
    real(dp) :: a(num_local * num_vecs), vecNorm(num_vecs), expect(num_vecs)

    OUT0("Starting TpetraMultiVector_get1dView!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec, 1.0_dp)
    a = vec%get1dView()

    TEST_FLOATING_ARRAY_EQUALITY(a, 1.0_dp, epsilon(1.0_dp))

    call vec%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_get1dView!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_get1dView)

  ! ----------------------------get1dViewNonConst----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_get1dViewNonConst)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2
    real(dp) :: a(num_local * num_vecs), vecNorm(num_vecs), expect(num_vecs)

    OUT0("Starting TpetraMultiVector_get1dViewNonConst!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec, 1.0_dp)
    a = vec%get1dViewNonConst()

    TEST_FLOATING_ARRAY_EQUALITY(a, 1.0_dp, epsilon(1.0_dp))

    call vec%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_get1dViewNonConst!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_get1dViewNonConst)

  ! --------------------------------sync_host--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_sync_host)
    type(TpetraMultiVector) :: Obj
    OUT0("Starting TpetraMultiVector_sync_host!")

    success = .false.

    !call Obj%create(); TEST_IERR()
    !call Obj%sync_host(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraMultiVector_sync_host: Test not yet implemented'

    OUT0("Finished TpetraMultiVector_sync_host!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_sync_host)

  ! -------------------------------sync_device-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_sync_device)
    type(TpetraMultiVector) :: Obj
    OUT0("Starting TpetraMultiVector_sync_device!")

    success = .false.

    !call Obj%create(); TEST_IERR()
    !call Obj%sync_device(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraMultiVector_sync_device: Test not yet implemented'

    OUT0("Finished TpetraMultiVector_sync_device!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_sync_device)

  ! ------------------------------need_sync_host------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_need_sync_host)
    type(TpetraMultiVector) :: Obj
    OUT0("Starting TpetraMultiVector_need_sync_host!")

    success = .false.

    !call Obj%create(); TEST_IERR()
    !fresult = Obj%need_sync_host(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraMultiVector_need_sync_host: Test not yet implemented'

    OUT0("Finished TpetraMultiVector_need_sync_host!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_need_sync_host)

  ! -----------------------------need_sync_device----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_need_sync_device)
    type(TpetraMultiVector) :: Obj
    OUT0("Starting TpetraMultiVector_need_sync_device!")
    success = .false.

    !call Obj%create(); TEST_IERR()
    !fresult = Obj%need_sync_device(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraMultiVector_need_sync_device: Test not yet implemented'

    OUT0("Finished TpetraMultiVector_need_sync_device!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_need_sync_device)

  ! ------------------------------modify_device------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_modify_device)
    type(TpetraMultiVector) :: Obj

    OUT0("Starting TpetraMultiVector_modify_device!")

    success = .false.

    !call Obj%create(); TEST_IERR()
    !fresult = Obj%modify_device(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraMultiVector_need_sync_device: Test not yet implemented'

    OUT0("Finished TpetraMultiVector_modify_device!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_modify_device)

  ! -------------------------------modify_host-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_modify_host)
    type(TpetraMultiVector) :: Obj

    OUT0("Starting TpetraMultiVector_modify_host!")

    success = .false.

    !call Obj%create(); TEST_IERR()
    !fresult = Obj%modify_host(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraMultiVector_need_sync_device: Test not yet implemented'

    OUT0("Finished TpetraMultiVector_modify_host!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_modify_host)

  ! ------------------------------getNumVectors------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_getNumVectors)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2

    OUT0("Starting TpetraMultiVector_getNumVectors!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    TEST_EQUALITY(num_vecs, vec%getNumVectors())
    call vec%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_getNumVectors!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_getNumVectors)

  ! ------------------------------getLocalLength------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_getLocalLength)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs =	2

    OUT0("Starting TpetraMultiVector_getLocalLength!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    TEST_EQUALITY(num_local, vec%getLocalLength())
    call vec%release();	TEST_IERR()

    OUT0("Finished TpetraMultiVector_getLocalLength!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_getLocalLength)

  ! -----------------------------getGlobalLength------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_getGlobalLength)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2

    OUT0("Starting TpetraMultiVector_getGlobalLength!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    TEST_EQUALITY(int(num_local * comm%getSize(), kind = size_type), vec%getGlobalLength())
    call vec%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_getGlobalLength!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_getGlobalLength)

  ! --------------------------------getStride--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_getStride)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2

    OUT0("Starting TpetraMultiVector_getStride!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    TEST_EQUALITY(int(num_local, kind = size_type), vec%getStride())
    call vec%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_getStride!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_getStride)

  ! -----------------------------isConstantStride----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_isConstantStride)
    type(TpetraMultiVector) :: vec
    integer, parameter :: num_local = 4
    integer(size_type), parameter :: num_vecs = 2

    OUT0("Starting TpetraMultiVector_isConstantStride!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec)
    TEST_ASSERT(vec%isConstantStride())
    call vec%release(); TEST_IERR()

    OUT0("Finished TpetraMultiVector_isConstantStride!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_isConstantStride)

  ! --------------------------------isSameSize-------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_isSameSize)
    type(TpetraMultiVector) :: vec1, vec2
    integer, parameter :: num_local = 4, bad_local = 6
    integer(size_type), parameter :: num_vecs = 2

    OUT0("Starting TpetraMultiVector_isSameSize!")

    call Tpetra_MV_Create(comm, num_local, num_vecs, vec1)
    call Tpetra_MV_Create(comm, num_local, num_vecs, vec2)
    TEST_ASSERT(vec1%isSameSize(vec2))

    call vec2%release(); TEST_IERR()
    call Tpetra_MV_Create(comm, bad_local, num_vecs, vec2)
    TEST_ASSERT(.not. vec1%isSameSize(vec2))

    call vec1%release(); TEST_IERR()
    call vec2%release(); TEST_IERR()
    OUT0("Finished TpetraMultiVector_isSameSize!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_isSameSize)

#if 0

  ! -------------------------------normWeighted------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_normWeighted)
    type(TpetraMultiVector) :: Obj
    type(TpetraMultiVector) :: weights
    real(scalar_type) :: norms(1)

    success = .false.

    !weights = xxx(); TEST_IERR()
    !norms = xxx(); TEST_IERR()
    !Obj = TpetraMultiVector(); TEST_IERR()
    !call Obj%normWeighted(weights, norms); TEST_IERR()
    !call weights%release(); TEST_IERR()
    !call norms%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'normWeighted: Test not yet implemented'

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_normWeighted)

  ! -----------------------removeEmptyProcessesInPlace------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_removeEmptyProcessesInPlace)
    type(TpetraMultiVector) :: Obj
    type(TpetraMap) :: newmap

    success = .false.

    !newmap = TpetraMap(); TEST_IERR()
    !Obj = TpetraMultiVector(); TEST_IERR()
    !call Obj%removeEmptyProcessesInPlace(newmap); TEST_IERR()
    !call newmap%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'removeEmptyProcessesInPlace: Test not yet implemented'

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_removeEmptyProcessesInPlace)

  ! ------------------------------setCopyOrView------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_setCopyOrView)
    type(TpetraMultiVector) :: Obj
    integer(kind(TeuchosDataAccess)) :: copyorview

    success = .false.

    copyorview = TeuchosCopy
    !Obj = TpetraMultiVector(); TEST_IERR()
    !call Obj%setCopyOrView(copyorview); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'setCopyOrView: Test not yet implemented'

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_setCopyOrView)

  ! ------------------------------getCopyOrView------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_getCopyOrView)
    type(TpetraMultiVector) :: Obj

    success = .false.

    !Obj = TpetraMultiVector(); TEST_IERR()
    !fresult = Obj%getCopyOrView(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'getCopyOrView: Test not yet implemented'

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_getCopyOrView)
#endif

end program test_TpetraMultiVector
