!Copyright 2017, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TpetraImportExport
#include "ForTrilinosTpetra_config.hpp"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  implicit none
  type(TeuchosComm) :: comm
  integer(global_size_type), parameter :: invalid=-1
  character(len=30), parameter :: FILENAME="test_tpetraimport_export.f90"

  SETUP_TEST()

#ifdef HAVE_MPI
  call comm%create(MPI_COMM_WORLD); CHECK_IERR()
#else
  call comm%create()
#endif

  ADD_SUBTEST_AND_RUN(TpetraImportExport_Basic)
  ADD_SUBTEST_AND_RUN(TpetraImportExport_GetNeighborsForward)

  !ADD_SUBTEST_AND_RUN(TpetraImport_isLocallyComplete)
  !ADD_SUBTEST_AND_RUN(TpetraImport_createRemoteOnlyImport)

  call comm%release()

  TEARDOWN_TEST()

contains

  ! -----------------------------------Basic --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_Basic)
    type(TpetraMap) :: src, tgt
    type(TpetraImport) :: importer
    integer(size_type), parameter :: ten=10, five=5
    integer(size_type) :: same, permute, remote, expected_sum, the_sum

    OUT0("Starting TpetraImportExport_Basic!")

    ! create Maps
    call src%create(invalid, ten, comm)
    call tgt%create(invalid, five, comm)

    ! create Import object
    call importer%create(src, tgt)

    same = importer%getNumSameIDs()
    permute = importer%getNumPermuteIDs()
    remote = importer%getNumRemoteIDs()
    the_sum = same + permute + remote
    expected_sum = tgt%getNodeNumElements()
    TEST_EQUALITY(the_sum, expected_sum)

    OUT0("Finished TpetraImportExport_Basic!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_Basic)

  ! ----------------------------- GetNeighborsForward ------------------------ !

  FORTRILINOS_UNIT_TEST(TpetraImportExport_GetNeighborsForward)
    type(TpetraMap) :: src, tgt
    type(TpetraImport) :: importer
    type(TpetraExport) :: exporter
    type(TpetraMultiVector) :: mv_mine, mv_with_neighbors
    type(TpetraMultiVector) :: mine_parent, neigh_parent
    real(scalar_type), allocatable :: a(:), val(:)
    real(scalar_type), parameter :: zero=0
    integer(size_type), parameter :: ten=10, five=5
    integer(size_type) :: num_images, my_image_id, num_local, num_vecs
    integer(size_type) :: tnum, j, n
    integer(local_ordinal_type) :: lclrow
    integer(global_ordinal_type), allocatable :: neighbors(:), cols(:)

    OUT0("Starting TpetraImportExport_GetNeighborsForward!")

    ! import with the importer to duplicate
    ! export with the exporter to add and reduce
    num_images = comm%getSize()
    my_image_id = comm%getRank()

    if (num_images < 2) return

    ! create a Map
    num_local = 1
    num_vecs = 5

    ! my neighbors: my_image_id-1, me, my_image_id+1
    if (my_image_id == 0 .or. my_image_id == num_images-1) then
      allocate(neighbors(2), a(2), val(2))
    else
      allocate(neighbors(3), a(3), val(3))
    end if

    if (my_image_id == 0) then
      neighbors(:) = [my_image_id+1, my_image_id+2]
    else if (my_image_id == num_images - 1) then
      neighbors(:) = [my_image_id, my_image_id+1]
    else
      neighbors(:) = [my_image_id, my_image_id+1, my_image_id+2]
    end if

    ! two maps: one has one entries per node, the other is the 1-D neighbors
    call src%create(invalid, num_local, comm)
    call tgt%create(invalid, neighbors, comm)

    do tnum = 1, 2

      ! for tnum=1, these are contiguously allocated multivectors
      ! for tnum=2, these are non-contiguous views of multivectors
      if (tnum == 1) then
        call mv_mine%create(src, num_vecs)
        call mv_with_neighbors%create(tgt, num_vecs)
      else
        allocate(cols(5))
        cols(:) = [1,7,4,5,6]
        call mine_parent%create(src, num_vecs+2)
        call neigh_parent%create(tgt, num_vecs+2)
        TEST_ASSERT((num_vecs == 5))
        call mv_mine%create()
        call mv_with_neighbors%create()
        mv_mine = mine_parent%subViewNonConst(cols)
        mv_with_neighbors = neigh_parent%subViewNonConst(cols)
      end if

      ! mv_mine = [my_image_id  my_image_id+num_images ... my_image_id+4*num_images]
      do j = 1, num_vecs
        lclrow = 1
        val(1) = real(my_image_id + j*num_images, kind=scalar_type)
        call mv_mine%replaceLocalValue(lclrow, j, val(1))
      end do

      ! create Import from src to tgt, Export from tgt to src, test them
      call importer%create(src, tgt)
      call exporter%create(tgt, src)

      ! importer testing
      TEST_ASSERT((src%isSameAs(importer%getSourceMap())))
      TEST_ASSERT((tgt%isSameAs(importer%getTargetMap())))

      n = merge(1, 0, my_image_id == 0)
      TEST_EQUALITY(importer%getNumSameIDs(), n)

      n = merge(0, 1, my_image_id == 0)
      TEST_EQUALITY(importer%getNumPermuteIDs(), n)

      n = merge(1, 2, (my_image_id == 0 .or. my_image_id == num_images - 1))
      TEST_EQUALITY(importer%getNumExportIDs(), n)
      TEST_EQUALITY(importer%getNumRemoteIDs(), n)

      ! exporter testing
      TEST_ASSERT(tgt%isSameAs(exporter%getSourceMap()))
      TEST_ASSERT(src%isSameAs(exporter%getTargetMap()))

      n = merge(1, 0, my_image_id == 0)
      TEST_EQUALITY(importer%getNumSameIDs(), n)

      n = merge(0, 1, my_image_id == 0)
      TEST_EQUALITY(exporter%getNumPermuteIDs(), n)

      ! import neighbors, test their proper arrival
      !                     [ 0    n     2n    3n    4n ]
      ! mv_with_neighbors = [...  ....  ....  ....  ....]
      !                     [n-1  2n-1  3n-1  4n-1  5n-1]
      call mv_with_neighbors%doImport(mv_mine, importer, TpetraREPLACE)
      do j = 1, num_vecs
        a = mv_with_neighbors%getData(j)
        if (my_image_id == 0) then
          val(1) = real(my_image_id+j*num_images,kind=scalar_type)
          val(2) = real(j*num_images+1, kind=scalar_type)
        else if (my_image_id == num_images-1) then
          val(1) = real(my_image_id+j*num_images-1, kind=scalar_type)
          val(2) = real(my_image_id+j*num_images, kind=scalar_type)
        else
          val(1) = real(my_image_id+j*num_images-1, kind=scalar_type)
          val(2) = real(my_image_id+j*num_images, kind=scalar_type)
          val(3) = real(my_image_id+j*num_images+1, kind=scalar_type)
        end if
        TEST_FLOATING_ARRAY_EQUALITY(a, val, epsilon(val(1))) ! me
        TEST_FLOATING_ARRAY_EQUALITY(a, val, epsilon(val(1))) ! neighbor
      end do

      ! export values, test
      call mv_mine%putScalar(zero)
      call mv_mine%doExport(mv_with_neighbors, exporter, TpetraADD)
      do j = 1, num_vecs
        a = mv_mine%getData(j)
        if (my_image_id == 0 .or. my_image_id == num_images-1) then
          ! contribution from me and one neighbor: double original value
          val(1) = real(2.0*(my_image_id+j*num_images), kind=scalar_type)
        else
          val(1) = real(3.0*(my_image_id+j*num_images), kind=scalar_type)
        end if
        TEST_FLOATING_EQUALITY(a(1), val(1), epsilon(val(1)))
      end do

    end do

    deallocate(neighbors, a, val)
    call src%release()
    call tgt%release()
    call mv_mine%release()
    call mv_with_neighbors%release()
    call mine_parent%release()
    call neigh_parent%release()
    call importer%release()
    call exporter%release()

    OUT0("Finished TpetraImportExport_GetNeighborsForward!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_GetNeighborsForward)

  ! --------------------------------- AbsMax --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_AbsMax)
    type(TpetraMap) :: smap, dmap
    type(TpetraImport) :: importer
    type(TpetraExport) :: exporter
    type(TpetraMultiVector) :: svec, dvec
    integer(size_type) :: num_images
    real(scalar_type) :: a(2)
    integer(local_ordinal_type) :: lclrow
    integer(global_ordinal_type) :: my_only_gid, cols(2)
    integer(size_type), parameter :: one=1
    integer(global_ordinal_type), parameter :: twog=2

    OUT0("Starting TpetraImportExport_AbsMax!")

    ! test ABSMAX CombineMode
    ! test with local and remote entries, as copyAndPermute() and
    ! unpackAndCombine() both need to be tested
    num_images = comm%getSize()

    if (num_images < 2) return;

    ! create a Map
    call smap%create(invalid, one, comm)

    lclrow = 1
    my_only_gid = smap%getGlobalElement(lclrow)
    cols(1) = my_only_gid
    cols(2) = mod(my_only_gid+1, num_images)
    call dmap%create(twog, cols, comm)

    call svec%create(smap, one)
    call svec%putScalar(-1.0_scalar_type)

    call dvec%create(dmap, one)
    call dvec%putScalar(-3.0_scalar_type)

    ! first item of dvec is local (w.r.t. srcVec), while the second is remote
    ! ergo, during the import:
    ! - the first will be over-written (by 1.0) from the source, while
    ! - the second will be "combined", i.e., abs(max(1.0,3.0)) = 3.0 from the dest
    call importer%create(smap, dmap)
    call dvec%doImport(svec, importer, TpetraABSMAX)

    a = dvec%get1dView()
    TEST_FLOATING_EQUALITY(a(1), -1.0_scalar_type, epsilon(a(1)))
    TEST_FLOATING_EQUALITY(a(2), 3.0_scalar_type, epsilon(a(1)))

    OUT0("Finished TpetraImportExport_AbsMax!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_AbsMax)

  ! ----------------------------isLocallyComplete----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImport_isLocallyComplete)
    type(TpetraImport) :: Obj
    OUT0("Starting TpetraImport_isLocallyComplete!")

    success = .false.

    !call Obj%create(); TEST_IERR()
    !fresult = Obj%isLocallyComplete(); TEST_IERR()

    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraImport_isLocallyComplete: Test not yet implemented'

    OUT0("Finished TpetraImport_isLocallyComplete!")

  END_FORTRILINOS_UNIT_TEST(TpetraImport_isLocallyComplete)

  ! --------------------------createRemoteOnlyImport-------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImport_createRemoteOnlyImport)
    type(TpetraImport) :: Obj
    type(TpetraMap) :: remotetarget
    OUT0("Starting TpetraImport_createRemoteOnlyImport!")

    success = .false.

    !call remotetarget%create(); TEST_IERR()
    !call Obj%create(); TEST_IERR()
    !fresult = Obj%createRemoteOnlyImport(remotetarget); TEST_IERR()

    !call remotetarget%release(); TEST_IERR()
    !call Obj%release(); TEST_IERR()

    write(*,*) 'TpetraImport_createRemoteOnlyImport: Test not yet implemented'

    OUT0("Finished TpetraImport_createRemoteOnlyImport!")

  END_FORTRILINOS_UNIT_TEST(TpetraImport_createRemoteOnlyImport)

end program test_TpetraImportExport
