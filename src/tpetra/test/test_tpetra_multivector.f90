program test_TpetraMultiVector
#include "ForTrilinosTpetra_config.hpp"
#include "FortranTestMacros.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  DECLARE_TEST_VARIABLES()
  type(TeuchosComm) :: comm
  integer(global_ordinal_type), parameter :: index_base = 1
  integer(global_size_type), parameter :: invalid = -1

  INITIALIZE_TEST()

#ifdef HAVE_MPI
  call comm%create(MPI_COMM_WORLD)
  CHECK_IERR()
#else
  call comm%create()
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
  CHECK_IERR()

  SHUTDOWN_TEST()

contains

  ! -----------------------------ZeroScaleUpdate------------------------------ !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_ZeroScaleUpdate)
    integer :: i
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: A, B, A2, C
    type(TeuchosArrayViewDoubleConst) :: AV
    integer(size_type), parameter :: num_vecs=2, num_local=2, LDA=2
    integer(local_ordinal_type) :: lclrow
    logical(bool_type) :: zeroout
    real(scalar_type), parameter :: zero=0., one=1., two=2., four=4., negone=-1.
    real(scalar_type) :: norms(num_vecs), zeros(num_vecs), values(6)
    integer(global_ordinal_type) :: gblrow, num_global

    zeros = zero

    call map%create(invalid, num_local, index_base, comm); TEST_IERR()

    ! values = {1, 1, 2, 2, 4, 4}
    ! values(1:4) = {1, 1, 2, 2} = [1 2]
    !                            = [1 2]
    ! values(3:)  = {2, 2, 4, 4} = [2 4]
    !                            = [2 4]
    ! a multivector A constructed from the first
    ! has values .5 of a multivector B constructed from the second
    ! then 2*A - B = 0
    ! we test both scale(), both update(), and norm()
    values = [one, one, two, two, four, four]

    ! TODO: Multivec create to take array, not ArrayView
    call AV%create(values(1:4)); TEST_IERR()
    call A%create(map, AV, LDA, num_vecs); TEST_IERR()
    call AV%release()

    call AV%create(values(3:)); TEST_IERR()
    call B%create(map, AV, LDA, num_vecs); TEST_IERR()
    call AV%release()

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
    call A2%create(A, Copy); TEST_IERR()
    call A2%scale(two)
    call A2%update(negone, B, one)
    call A2%norm1(norms)
    TEST_COMPARE_FLOATING_ARRAYS(norms, zeros, epsilon(zero))
    call A2%release()

    ! set A2 = A
    ! check that it equals B: scale, subtraction in situ
    call A2%create(A, Copy); TEST_IERR()
    call A2%update(negone, B, two)
    call A2%norm1(norms)
    TEST_COMPARE_FLOATING_ARRAYS(norms, zeros, epsilon(zero))
    call A2%release()

    ! set C random
    ! set it to zero by combination with A,B
    zeroout = .false.
    call C%create(map, num_vecs, zeroout); TEST_IERR()
    call C%randomize()
    call C%update(negone, B, two, A, zero)
    call C%norm1(norms)
    TEST_COMPARE_FLOATING_ARRAYS(norms, zeros, epsilon(zero))
    call C%release()

    ! set C random
    ! scale it ex-situ
    ! check that it equals B: subtraction in situ
    call C%create(map, num_vecs, zeroout); TEST_IERR()
    call C%scale(two, A)
    call C%update(one, B, negone)
    call C%norm1(norms)
    TEST_COMPARE_FLOATING_ARRAYS(norms, zeros, epsilon(zero))
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
    type(TeuchosArrayViewDoubleConst) :: AV
    integer(size_type), parameter :: num_vecs=3, num_local=2, LDA=2
    real(scalar_type) :: values(num_vecs*num_local), answer(num_vecs), norms(num_vecs)

    OUT0("Starting CountNormInf")

    ! create a Map
    call map%create(invalid, num_local, index_base, comm); TEST_IERR()
    ! values = {0, 0, 1, 1, 2, 2} = [0 1 2]
    !                               [0 1 2]
    ! normInf(values) = [0 1 2]
    ! over all procs, this is [0 1 2]
    values(:) = [0., 0., 1., 1., 2., 2.]
    ! TODO: MultiVector takes array instead of ArrayView
    call AV%create(values); TEST_IERR()
    call Vec%create(map, AV, LDA, num_vecs); TEST_IERR()
    call AV%release()
    answer(:) = [0., 1., 2.]

    ! do the dots
    call Vec%normInf(norms)

    call Vec%release()
    call map%release()

    ! check the answers
    TEST_COMPARE_FLOATING_ARRAYS(norms, answer, epsilon(answer(1)))

    OUT0("Finished CountNormInf!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_CountNormInf)

  ! --------------------------------Norm2------------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Norm2)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    real(scalar_type), parameter :: zero=0.
    integer(size_type), parameter :: num_vecs=7, num_local=13
    real(scalar_type) :: norms_rand(num_vecs), norms_zero(num_vecs)

    OUT0("Starting Norm2")

    ! create a Map
    call map%create(invalid, num_local, index_base, comm); TEST_IERR()

    call Vec%create(map, num_vecs); TEST_IERR()
    call Vec%randomize(); TEST_IERR()

    ! Take the norms, they should not be zero
    call Vec%norm2(norms_rand)

    ! Zero the vector
    call Vec%putScalar(zero)

    ! Take the norms, they should be zero
    call Vec%norm2(norms_zero)

    ! Check the answers
    TEST_ARRAY_INEQUALITY(norms_rand, zero, epsilon(zero))
    TEST_ARRAY_EQUALITY(norms_zero, zero, epsilon(zero))

    call Vec%release()
    call map%release()

    OUT0("Finished Norm2!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Norm2)

  ! -----------------------------Reciprocal----------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Reciprocal)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: A, B
    integer(size_type), parameter :: num_vecs=2, num_local=10
    real(scalar_type) :: dots(num_vecs)
    real(scalar_type), parameter :: one=1., five=5.

    OUT0("Starting Reciprocal")

    ! create a Map
    call map%create(invalid, num_local, index_base, comm); TEST_IERR()

    call A%create(map, num_vecs); TEST_IERR()
    call B%create(map, num_vecs); TEST_IERR()
    call A%putScalar(five); TEST_IERR()
    call A%reciprocal(B); TEST_IERR()

    ! Take the dots, they should one
    call A%dot(B, dots)

    ! Check the answers
    TEST_ARRAY_EQUALITY(dots, one, epsilon(one))

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
    call map1%create(num_global, index_base, comm); TEST_IERR()

    num_global = 5 * comm%getSize()
    call map2%create(num_global, index_base, comm); TEST_IERR()

    call Obj%create(map1, num_vecs); TEST_IERR()

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
    integer(size_type), parameter :: num_vecs=2, num_local=10
    real(scalar_type), parameter :: zero=0., one=1., negone=-1.
    real(scalar_type) :: norms(num_vecs)

    OUT0("Starting Abs")
    call map%create(invalid, num_local, index_base, comm); TEST_IERR()

    call A%create(map, num_vecs); TEST_IERR()
    call B%create(map, num_vecs); TEST_IERR()
    call A%putScalar(negone); TEST_IERR()
    call A%abs(B)

    !   set A2 = A
    !   scale it by 2 in situ
    !   check that it equals B: subtraction in situ
    call A2%create(A, Copy); TEST_IERR()
    call A2%update(one, B, one)
    call A2%norm1(norms)
    TEST_ARRAY_EQUALITY(norms, zero, epsilon(zero))

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
    type(string) :: fresult
    integer(size_type), parameter :: num_vecs=2, num_local=10
    integer(global_ordinal_type) :: num_global, index_base
    call map%create(invalid, num_local, index_base, comm); TEST_IERR()
    call Vec%create(map, num_vecs)
    fresult = Vec%description(); TEST_IERR()
    call Vec%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Description)

  ! --------------------------------MeanValue--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_MeanValue)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    type(TeuchosArrayViewDoubleConst) :: AV
    type(TeuchosArrayViewDouble) :: means_av
    integer(size_type), parameter :: num_vecs=2, num_local=2, LDA=2
    real(scalar_type) :: values(4), means(num_vecs), answer(num_vecs)

    OUT0("Starting MeanValue")

    call map%create(invalid, num_local, index_base, comm); TEST_IERR()

    ! values = {2, 6, 3, 1} = [2 3]
    !                         [6 1]
    values(:) = [2., 6., 3., 1.]

    ! TODO: Multivec create to take array, not ArrayView
    call AV%create(values); TEST_IERR()
    call Vec%create(map, AV, LDA, num_vecs); TEST_IERR()
    call AV%release()

    !call means%create(); TEST_IERR()
    call means_av%create(means); TEST_IERR()
    call Vec%meanValue(means_av); TEST_IERR()
    call means_av%release(); TEST_IERR()

    answer = [4., 2.]
    TEST_COMPARE_FLOATING_ARRAYS(means, answer, epsilon(means(1)))

    call Vec%release(); TEST_IERR()
    call map%release(); TEST_IERR()

    OUT0("Finished MeanValue")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_MeanValue)

  ! ---------------------------------Multiply--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Multiply)
    type(TpetraMap) :: map2, map3
    type(TpetraMultiVector) :: Vec2x2, Vec3x2
    real(scalar_type), parameter :: zero=0., one=1.
    integer(size_type), parameter :: n2=2, n3=3
    integer(int_type) :: num_images
    integer(global_size_type) :: num_global
    real(scalar_type) :: check(9)

    OUT0("Starting Multiply")

    num_images = comm%getSize()
    check = 3 * num_images

    num_global = n2 * comm%getSize()
    call map2%create(num_global, index_base, comm, LocallyReplicated); TEST_IERR()
    call map3%create(invalid, n3, index_base, comm); TEST_IERR()

    call Vec2x2%create(map2, n2); TEST_IERR()
    call Vec3x2%create(map3, n2); TEST_IERR()
    call Vec3x2%putScalar(one); TEST_IERR()

    call Vec2x2%multiply(CONJ_TRANS, NO_TRANS, one, Vec3x2, Vec3x2, zero)
    TEST_IERR()

    !a = fortran array view of the data
    !n = size of a
    !TEST_COMPARE_FLOATING_ARRAYS(a, check(1:n), epsilon(a(1)))

    OUT0("Finished Multiply!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Multiply)

  ! --------------------------------Basic------------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Basic)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    integer(size_type), parameter :: num_vecs=12, num_local=2
    call map%create(invalid, num_local, index_base, comm); TEST_IERR()
    call Vec%create(map, num_vecs); TEST_IERR()

    TEST_EQUALITY(Vec%getNumVectors(), num_vecs)
    TEST_EQUALITY(Vec%getLocalLength(), num_local)
    TEST_EQUALITY(Vec%getGlobalLength(), num_local*comm%getSize())
    TEST_EQUALITY(Vec%getStride(), num_local)
    TEST_ASSERT(Vec%isConstantStride())

    call Vec%release(); TEST_IERR()
    call map%release(); TEST_IERR()
  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Basic)

  ! ----------------------------------Reduce---------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Reduce)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    integer(global_size_type) :: num_global
    integer(size_type), parameter :: num_vecs=1, num_local=2
    real(scalar_type), parameter :: two=2.
    OUT0("Starting Reduce")
    num_global = num_local * comm%getSize()
    call map%create(num_global, index_base, comm, LocallyReplicated); TEST_IERR()
    call Vec%create(map, num_vecs); TEST_IERR()
    call Vec%putScalar(two); TEST_IERR()
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
    integer(size_type), parameter :: num_vecs=2, num_local=4
    integer(local_ordinal_type) :: lclrow
    real(scalar_type) :: value, expected, dots(num_vecs)
    real(scalar_type), parameter :: one=1.
    integer(global_ordinal_type) :: gblrow, num_global
    OUT0("Starting replaceGlobalValue")
    num_global = num_local * comm%getSize()
    call map%create(num_global, index_base, comm); TEST_IERR()

    call Vec%create(map, num_vecs); TEST_IERR()
    call OneV%create(map, num_vecs); TEST_IERR()
    call OneV%putScalar(one); TEST_IERR()

    ! TODO: Columns should be 1 based
    do lclrow = 1, num_local
      gblrow = map%getGlobalElement(lclrow)
      value = real(gblrow, kind=scalar_type)
      do col = 0, num_vecs-1
        call Vec%replaceGlobalValue(gblrow, col, value); TEST_IERR()
      end do
    end do

    call Vec%dot(OneV, dots)
    expected = real(num_global * (num_global + 1) / 2., kind=scalar_type)
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
    integer(local_ordinal_type) :: lclrow
    integer(size_type), parameter :: num_vecs=2, num_local=4
    real(scalar_type) :: value, dots(num_vecs), expected
    real(scalar_type), parameter :: one=1.
    OUT0("Starting ReplaceLocalValue")

    call map%create(invalid, num_local, index_base, comm); TEST_IERR()
    call Vec%create(map, num_vecs); TEST_IERR()
    call OneV%create(map, num_vecs); TEST_IERR()
    call OneV%putScalar(one)

    ! TODO: column indices are wrong!  The are 0 based, but there is something
    ! wrong with the stride, uncomment the assertions below to see.
    do lclrow = 1, num_local
      value = real(lclrow, kind=scalar_type)
      do col = 0, num_vecs-1
        call Vec%replaceLocalValue(lclrow, col, value); TEST_IERR()
      end do
    end do

    call Vec%dot(OneV, dots)
    expected = real(comm%getSize() * (num_local * (num_local + 1)) / 2., kind=scalar_type)
    !TODO: print*, expected
    !TODO: print*, dots
    !TODO: TEST_FLOATING_EQUALITY(dots(1), dots(2), epsilon(dots(2)))
    !TODO: TEST_FLOATING_EQUALITY(expected, dots(1), epsilon(dots(2)))

    call Vec%release(); TEST_IERR()
    call OneV%release(); TEST_IERR()
    call map%release()

    OUT0("Finished ReplaceLocalValue!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_ReplaceLocalValue)

  ! --------------------------------Get1dCopy--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_Get1dCopy)
    type(TpetraMap) :: map
    type(TpetraMultiVector) :: Vec
    type(TeuchosArrayViewDouble) :: AV
    integer(size_type) :: lda
    real(scalar_type), allocatable :: a(:)
    integer(size_type), parameter :: num_vecs=2, num_local=4
    real(scalar_type), parameter :: one=1.
    OUT0("Starting Get1dCopy")
    call map%create(invalid, num_local, index_base, comm); TEST_IERR()
    call Vec%create(map, num_vecs); TEST_IERR()
    call Vec%putScalar(one); TEST_IERR()
    allocate(a(num_vecs*num_local*comm%getSize()))
    a = 0.
    call AV%create(a); TEST_IERR()
    lda = num_local*comm%getSize()
    call Vec%Get1dCopy(AV, lda); TEST_IERR() ! TODO: Get1dCopy to take fortran array
    call AV%release(); TEST_IERR()
    call Vec%release(); TEST_IERR()

    TEST_ARRAY_EQUALITY(a, one, epsilon(one))
    OUT0("Finished Get1dCopy!")

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_Get1dCopy)


  ! ----------------------------offsetViewNonConst---------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_offsetViewNonConst)
    type(TpetraMultiVector) :: Obj
    type(TpetraMap) :: submap
    integer(size_type) :: offset

    success = .false.

    !call submap%create(); TEST_IERR()
    offset = 0
    !call Obj%create(); TEST_IERR()
    !fresult = Obj%offsetViewNonConst(submap, offset); TEST_IERR()
    !call submap%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'offsetViewNonConst: Test not yet implemented'

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_offsetViewNonConst)

  ! -------------------------------normWeighted------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_normWeighted)
    type(TpetraMultiVector) :: Obj
    type(TpetraMultiVector) :: weights
    type(TeuchosArrayViewDouble) :: norms

    success = .false.

    !call weights%create(); TEST_IERR()
    !call norms%create(); TEST_IERR()
    !call Obj%create(); TEST_IERR()
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

    !call newmap%create(); TEST_IERR()
    !call Obj%create(); TEST_IERR()
    !call Obj%removeEmptyProcessesInPlace(newmap); TEST_IERR()
    !call newmap%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'removeEmptyProcessesInPlace: Test not yet implemented'

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_removeEmptyProcessesInPlace)

  ! ------------------------------setCopyOrView------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_setCopyOrView)
    type(TpetraMultiVector) :: Obj
    integer(kind(DataAccess)) :: copyorview

    success = .false.

    copyorview = Copy
    !call Obj%create(); TEST_IERR()
    !call Obj%setCopyOrView(copyorview); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'setCopyOrView: Test not yet implemented'

  END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_setCopyOrView)

  ! ------------------------------getCopyOrView------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraMultiVector_getCopyOrView)
    type(TpetraMultiVector) :: Obj

    success = .false.

    !call Obj%create(); TEST_IERR()
    !fresult = Obj%getCopyOrView(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'getCopyOrView: Test not yet implemented'

    END_FORTRILINOS_UNIT_TEST(TpetraMultiVector_getCopyOrView)


end program test_TpetraMultiVector
