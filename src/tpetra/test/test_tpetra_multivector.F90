! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program test_TpetraMultiVector
#include "ForTrilinosTpetra_config.hpp"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  implicit none
  type(TeuchosComm) :: comm
  character(len=256), parameter :: FILENAME="test_tpetra_multivector.f90"

  SETUP_TEST()

#ifdef HAVE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm(); FORTRILINOS_CHECK_IERR()
#endif

  ADD_SUBTEST_AND_RUN(TpetraMultiVector_ZeroScaleUpdate)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_CountNormInf)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_Norm2)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_ReplaceMap)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_Reciprocal)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_Abs)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_Description)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_MeanValue)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_Multiply)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_Basic)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_Reduce)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_ReplaceGlobalValue)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_ReplaceLocalValue)
  ADD_SUBTEST_AND_RUN(TpetraMultiVector_Get1dCopy)

  ! TODO: The following tests have only skeletons
  !ADD_SUBTEST_AND_RUN(TpetraMultiVector_offsetViewNonConst)

  ! The following methods are deprecated
  !ADD_SUBTEST_AND_RUN(TpetraMultiVector_normWeighted)

  ! The following methods are listed as "may change or disappear at any time"
  !ADD_SUBTEST_AND_RUN(TpetraMultiVector_setCopyOrView)
  !ADD_SUBTEST_AND_RUN(TpetraMultiVector_getCopyOrView)
  !ADD_SUBTEST_AND_RUN(TpetraMultiVector_removeEmptyProcessesInPlace)

  call comm%release()

  TEARDOWN_TEST()

contains

  ! -----------------------------ZeroScaleUpdate------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_ZeroScaleUpdate)
    integer :: i
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: A, B, A2, C
    integer(size_type), parameter :: num_vecs=2, LDA=2
    integer, parameter :: num_local=2
    integer :: lclrow
    logical :: zeroout
    real(scalar_type) :: norms(num_vecs), zeros(num_vecs), values(6)
    integer(global_ordinal_type) :: gblrow, num_global

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
    call A2%release()

    ! set A2 = A
    ! check that it equals B: scale, subtraction in situ
    A2 = TpetraMultiVector(A, TeuchosCopy); TEST_IERR()
    call A2%update(-1.d0, B, 2.d0)
    call A2%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, zeros, epsilon(0.d0))
    call A2%release()

    ! set C random
    ! set it to zero by combination with A,B
    zeroout = .false.
    C = TpetraMultiVector(map, num_vecs, zeroout); TEST_IERR()
    call C%randomize()
    call C%update(-1.d0, B, 2.d0, A, 0.d0)
    call C%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, zeros, epsilon(0.d0))
    call C%release()

    ! set C random
    ! scale it ex-situ
    ! check that it equals B: subtraction in situ
    C = TpetraMultiVector(map, num_vecs, zeroout); TEST_IERR()
    call C%scale(2.d0, A)
    call C%update(1.0d0, B, -1.d0)
    call C%norm1(norms)
    TEST_FLOATING_ARRAY_EQUALITY(norms, zeros, epsilon(0.d0))
    call C%release()

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

    call Vec%release()
    call map%release()

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

    call Vec%release()
    call map%release()

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

    call A%release()
    call B%release()
    call map%release()

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
    integer :: i
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

    call A2%release()
    call B%release()
    call A%release()
    call map%release()

    OUT0("Finished Abs")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Abs)

! --------------------------------description--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Description)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    character(kind=C_CHAR, len=:), allocatable :: description
    integer(size_type), parameter :: num_vecs=2
    integer, parameter :: num_local=10
    integer(global_ordinal_type) :: num_global
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
    integer(global_size_type) :: num_global
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

  ! --------------------------------Basic------------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Basic)
    type(TpetraMap) :: map
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
    integer :: i
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
    call map%release()

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


#if 0
  ! ----------------------------offsetViewNonConst---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_offsetViewNonConst)
    type(TpetraMultiVector) :: Obj
    type(TpetraMap) :: submap
    integer(size_type) :: offset

    success = .false.

    !submap = TpetraMap(); TEST_IERR()
    offset = 0
    !Obj = TpetraMultiVector(); TEST_IERR()
    !fresult = Obj%offsetViewNonConst(submap, offset); TEST_IERR()
    !call submap%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'offsetViewNonConst: Test not yet implemented'

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_offsetViewNonConst)

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
