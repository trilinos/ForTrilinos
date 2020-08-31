!Copyright 2017-2018, UT-Battelle, LLC
!
!SPDX-License-Identifier: BSD-3-Clause
!License-Filename: LICENSE
program test_TpetraImportExport
#include "ForTrilinos_config.h"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra
  use test_tpetra_import_export_helper

  implicit none
  type(TeuchosComm) :: comm
  character(len=30), parameter :: FILENAME="test_tpetra_import_export.F90"

  SETUP_TEST()

#if FORTRILINOS_USE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm(); FORTRILINOS_CHECK_IERR()
#endif

  ADD_SUBTEST_AND_RUN(TpetraImportExport_Basic)
  ADD_SUBTEST_AND_RUN(TpetraImportExport_GetNeighborsForward)
  ADD_SUBTEST_AND_RUN(TpetraImportExport_AbsMax)

  ADD_SUBTEST_AND_RUN(TpetraImportExport_getNumSameIDs)
  ADD_SUBTEST_AND_RUN(TpetraImportExport_getNumRemoteIDs)
  ADD_SUBTEST_AND_RUN(TpetraImportExport_getNumPermuteIDs)
  ADD_SUBTEST_AND_RUN(TpetraImportExport_getNumExportIDs)
  ADD_SUBTEST_AND_RUN(TpetraImportExport_getSourceMap)
  ADD_SUBTEST_AND_RUN(TpetraImportExport_getTargetMap)

  ! Methods will not work with 1 rank so skip
  if(comm%getSize() /= 1) then
    ADD_SUBTEST_AND_RUN(TpetraImport_createRemoteOnlyImport)
  end if

  ! Obsolete test
  !ADD_SUBTEST_AND_RUN(TpetraImport_isLocallyComplete)

  call comm%release();   TEST_IERR()

  TEARDOWN_TEST()

contains

  ! -----------------------------------Basic --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_Basic)
    type(TpetraMap) :: tgt
    type(TpetraImport) :: importer
    type(TpetraExport) :: exporter
    integer(size_type) :: same, permute, remote, expected_sum, the_sum
    integer :: num_images

    OUT0("Starting TpetraImportExport_Basic!")

    num_images = comm%getSize()

    ! create Import object
    call Tpetra_IE_MakeImport(comm, 10 * num_images, 5 * num_images, importer)

    tgt = importer%getTargetMap()

    same = importer%getNumSameIDs()
    permute = importer%getNumPermuteIDs()
    remote = importer%getNumRemoteIDs()
    the_sum = same + permute + remote
    expected_sum = tgt%getNodeNumElements()
    TEST_EQUALITY(the_sum, expected_sum)

    ! Create Export and perform similar examination
    call tgt%release(); TEST_IERR()
    call Tpetra_IE_MakeExport(comm, 10 * num_images, 5 * num_images, exporter)

    tgt = exporter%getTargetMap()

    same = exporter%getNumSameIDs()
    permute = exporter%getNumPermuteIDs()
    remote = exporter%getNumRemoteIDs()
    the_sum = same + permute + remote
    expected_sum = tgt%getNodeNumElements()
    TEST_EQUALITY(the_sum, expected_sum)

    call tgt%release(); TEST_IERR()
    call importer%release(); TEST_IERR()
    call exporter%release(); TEST_IERR()

    OUT0("Finished TpetraImportExport_Basic!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_Basic)

  ! ----------------------------- GetNeighborsForward ------------------------ !

  FORTRILINOS_UNIT_TEST(TpetraImportExport_GetNeighborsForward)
    type(TpetraMap) :: src, tgt, tmp
    type(TpetraImport) :: importer
    type(TpetraExport) :: exporter
    type(TpetraMultiVector) :: mv_mine, mv_with_neighbors
    type(TpetraMultiVector) :: mine_parent, neigh_parent
    real(scalar_type), allocatable :: val(:)
    real(scalar_type), pointer :: a(:)
    integer(size_type) :: num_images, my_image_id, num_vecs
    integer(size_type) :: tnum, j, n
    integer :: lclrow, num_local
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
    src = TpetraMap(TPETRA_GLOBAL_INVALID, num_local, comm)
    tgt = TpetraMap(TPETRA_GLOBAL_INVALID, neighbors, comm)

    do tnum = 1, 2

      ! for tnum=1, these are contiguously allocated multivectors
      ! for tnum=2, these are non-contiguous views of multivectors
      if (tnum == 1) then
        mv_mine = TpetraMultiVector(src, num_vecs); TEST_IERR()
        mv_with_neighbors = TpetraMultiVector(tgt, num_vecs); TEST_IERR()
      else
        allocate(cols(5))
        cols(:) = [1,7,4,5,6]
        mine_parent = TpetraMultiVector(src, num_vecs+2); TEST_IERR()
        neigh_parent = TpetraMultiVector(tgt, num_vecs+2); TEST_IERR()
        TEST_ASSERT((num_vecs == 5)); TEST_IERR()
        mv_mine = mine_parent%subViewNonConst(cols); TEST_IERR()
        mv_with_neighbors = neigh_parent%subViewNonConst(cols); TEST_IERR()
      end if

      ! mv_mine = [my_image_id  my_image_id+num_images ... my_image_id+4*num_images]
      do j = 1, num_vecs
        lclrow = 1
        val(1) = real(my_image_id + j*num_images, kind=scalar_type)
        call mv_mine%replaceLocalValue(lclrow, j, val(1))
      end do

      ! create Import from src to tgt, Export from tgt to src, test them
      importer = TpetraImport(src, tgt)
      exporter = TpetraExport(tgt, src)

      ! importer testing
      tmp = importer%getSourceMap()
      TEST_ASSERT(src%isSameAs(tmp))
      call tmp%release();   TEST_IERR()
      tmp = importer%getTargetMap()
      TEST_ASSERT(tgt%isSameAs(tmp))
      call tmp%release();   TEST_IERR()

      n = merge(1, 0, my_image_id == 0)
      TEST_EQUALITY(importer%getNumSameIDs(), n)

      n = merge(0, 1, my_image_id == 0)
      TEST_EQUALITY(importer%getNumPermuteIDs(), n)

      n = merge(1, 2, (my_image_id == 0 .or. my_image_id == num_images - 1))
      TEST_EQUALITY(importer%getNumExportIDs(), n)
      TEST_EQUALITY(importer%getNumRemoteIDs(), n)

      ! exporter testing
      tmp = exporter%getSourceMap()
      TEST_ASSERT(tgt%isSameAs(tmp))
      call tmp%release();   TEST_IERR()
      tmp = exporter%getTargetMap()
      TEST_ASSERT(src%isSameAs(tmp))
      call tmp%release();   TEST_IERR()

      n = merge(1, 0, my_image_id == 0)
      TEST_EQUALITY(exporter%getNumSameIDs(), n)

      n = merge(0, 1, my_image_id == 0)
      TEST_EQUALITY(exporter%getNumPermuteIDs(), n)

      ! import neighbors, test their proper arrival
      !                     [ 0    n     2n    3n    4n ]
      ! mv_with_neighbors = [...  ....  ....  ....  ....]
      !                     [n-1  2n-1  3n-1  4n-1  5n-1]
      call mv_with_neighbors%doImport(mv_mine, importer, TpetraREPLACE)
      do j = 1, num_vecs
        a => mv_with_neighbors%getData(j)
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
      call mv_mine%putScalar(0d0)
      call mv_mine%doExport(mv_with_neighbors, exporter, TpetraADD)
      do j = 1, num_vecs
        a => mv_mine%getData(j)
        if (my_image_id == 0 .or. my_image_id == num_images-1) then
          ! contribution from me and one neighbor: double original value
          val(1) = real(2*(my_image_id+j*num_images), kind=scalar_type)
        else
          val(1) = real(3*(my_image_id+j*num_images), kind=scalar_type)
        end if
        TEST_FLOATING_EQUALITY(a(1), val(1), epsilon(val(1)))
      end do

      call mv_mine%release();   TEST_IERR()
      call mv_with_neighbors%release();   TEST_IERR()
      call mine_parent%release();   TEST_IERR()
      call neigh_parent%release();   TEST_IERR()
      call importer%release();   TEST_IERR()
      call exporter%release();   TEST_IERR()

    end do

    deallocate(neighbors, val)
    call src%release();   TEST_IERR()
    call tgt%release();   TEST_IERR()

    OUT0("Finished TpetraImportExport_GetNeighborsForward!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_GetNeighborsForward)

  ! --------------------------------- AbsMax --------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_AbsMax)
    type(TpetraMap) :: smap, dmap
    type(TpetraImport) :: importer
    type(TpetraMultiVector) :: svec, dvec
    integer(size_type) :: num_images
    real(scalar_type), pointer :: a(:)
    integer :: k
    integer(global_ordinal_type) :: num_global
    integer(global_ordinal_type), allocatable :: cols(:)
    integer(size_type), parameter :: one=1
    integer(global_ordinal_type), parameter :: twog=2

    OUT0("Starting TpetraImportExport_AbsMax!")

    ! test ABSMAX CombineMode
    ! test with local and remote entries, as copyAndPermute() and
    ! unpackAndCombine() both need to be tested
    num_images = comm%getSize()

    if (num_images < 2) return;
    allocate(cols(2))
    num_global = int(twog * num_images, kind = global_ordinal_type)

    ! create a Map
    smap = TpetraMap(TPETRA_GLOBAL_INVALID, 1, comm)

    !Use noncontiguous map setup from isContiguous() in test_tpetra_map
    do k = 1, 2
        cols(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do

    !Offset by 1 to make second element of dvec remote w.r.t svec
    cols(1) = int(comm%getRank()+1, kind=global_ordinal_type)
    dmap = TpetraMap(TPETRA_GLOBAL_INVALID, cols, comm)

    svec = TpetraMultiVector(smap, one)
    call svec%putScalar(-1.0_scalar_type)

    dvec = TpetraMultiVector(dmap, one)
    call dvec%putScalar(3.0_scalar_type)

    ! first item of dvec is local (w.r.t. svec), while the second is remote
    ! ergo, during the import:
    ! - the first will be over-written (by -1.0) from the source, while
    ! - the second will be "combined", i.e., abs(max(-1.0,3.0)) = 3.0 from the dest
    importer = TpetraImport(smap, dmap)
    call dvec%doImport(svec, importer, TpetraABSMAX)

    a => dvec%get1dView()
    TEST_FLOATING_EQUALITY(a(1), -1.0_scalar_type, epsilon(a(1)))
    TEST_FLOATING_EQUALITY(a(2), 3.0_scalar_type, epsilon(a(2)))

    deallocate(cols)
    call importer%release(); TEST_IERR()
    call svec%release(); TEST_IERR()
    call dvec%release(); TEST_IERR()
    call smap%release(); TEST_IERR()
    call dmap%release(); TEST_IERR()

    OUT0("Finished TpetraImportExport_AbsMax!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_AbsMax)

 ! ------------------------------getNumSameIDs------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_getNumSameIDs)
    type(TpetraImport) :: imp
    type(TpetraExport) :: exp
    type(TpetraMap) :: tgt
    integer :: numSame, numRemote, numPermute, numNodes
    OUT0("Starting TpetraImportExport_getNumSameIDs!")

    ! Use property that numSame + numRemote + numPermute = # elements on target map

    ! Import test
    call Tpetra_IE_MakeImport(comm, 5 * comm%getSize(), 2 * comm%getSize(), imp)
    tgt = imp%getTargetMap()
    numSame = imp%getNumSameIDs()
    numRemote = imp%getNumRemoteIDs()
    numPermute = imp%getNumPermuteIDs()
    numNodes = tgt%getNodeNumElements()
    TEST_EQUALITY(numSame, numNodes - numPermute - numRemote)

    ! Export test
    call tgt%release(); TEST_IERR()
    call Tpetra_IE_MakeExport(comm, 5 * comm%getSize(), 2 * comm%getSize(), exp)
    tgt = exp%getTargetMap()
    numSame = exp%getNumSameIDs()
    numRemote =	exp%getNumRemoteIDs()
    numPermute = exp%getNumPermuteIDs()
    numNodes = tgt%getNodeNumElements()
    TEST_EQUALITY(numSame, numNodes - numPermute - numRemote)

    call imp%release(); TEST_IERR()
    call exp%release(); TEST_IERR()
    call tgt%release(); TEST_IERR()

    OUT0("Finished TpetraImportExport_getNumSameIDs!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_getNumSameIDs)

 ! ------------------------------getNumRemoteIDs------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_getNumRemoteIDs)
    type(TpetraImport) :: imp
    type(TpetraExport) :: exp
    type(TpetraMap) :: tgt
    integer :: numSame, numRemote, numPermute, numNodes
    OUT0("Starting TpetraImportExport_getNumRemoteIDs!")

    ! Use property that numSame + numRemote + numPermute = # elements on target map

    ! Import test
    call Tpetra_IE_MakeImport(comm, 5 * comm%getSize(), 2 * comm%getSize(), imp)
    tgt = imp%getTargetMap()
    numSame = imp%getNumSameIDs()
    numRemote = imp%getNumRemoteIDs()
    numPermute = imp%getNumPermuteIDs()
    numNodes = tgt%getNodeNumElements()
    TEST_EQUALITY(numRemote, numNodes - numPermute - numSame)

    ! Export test
    call Tpetra_IE_MakeExport(comm, 5 * comm%getSize(), 2 * comm%getSize(), exp)
    call tgt%release(); TEST_IERR()
    numSame = exp%getNumSameIDs()
    numRemote = exp%getNumRemoteIDs()
    numPermute = exp%getNumPermuteIDs()
    tgt = exp%getTargetMap()
    numNodes = tgt%getNodeNumElements()
    TEST_EQUALITY(numRemote, numNodes - numPermute - numSame)

    call imp%release(); TEST_IERR()
    call exp%release(); TEST_IERR()
    call tgt%release(); TEST_IERR()

    OUT0("Finished TpetraImportExport_getNumRemoteIDs!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_getNumRemoteIDs)

 ! ------------------------------getNumPermuteIDs------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_getNumPermuteIDs)
    type(TpetraImport) :: imp
    type(TpetraExport) :: exp
    type(TpetraMap) :: tgt
    integer :: numSame, numRemote, numPermute, numNodes
    OUT0("Starting TpetraImportExport_getNumPermuteIDs!")

    ! Use property that numSame + numRemote + numPermute = # elements on target map

    ! Import test
    call Tpetra_IE_MakeImport(comm, 5 * comm%getSize(), 2 * comm%getSize(), imp)
    tgt = imp%getTargetMap()
    numSame = imp%getNumSameIDs()
    numRemote = imp%getNumRemoteIDs()
    numPermute = imp%getNumPermuteIDs()
    numNodes = tgt%getNodeNumElements()
    TEST_EQUALITY(numPermute, numNodes - numSame - numRemote)

    ! Export test
    call tgt%release(); TEST_IERR()
    call Tpetra_IE_MakeExport(comm, 5 * comm%getSize(), 2 * comm%getSize(), exp)
    numSame = exp%getNumSameIDs()
    numRemote = exp%getNumRemoteIDs()
    numPermute = exp%getNumPermuteIDs()
    tgt = exp%getTargetMap()
    numNodes = tgt%getNodeNumElements()
    TEST_EQUALITY(numPermute, numNodes - numSame - numRemote)

    call imp%release(); TEST_IERR()
    call exp%release(); TEST_IERR()
    call tgt%release(); TEST_IERR()

    OUT0("Finished TpetraImportExport_getNumPermuteIDs!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_getNumPermuteIDs)

 ! ------------------------------getNumExportIDs------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_getNumExportIDs)
    type(TpetraImport) :: imp
    type(TpetraExport) :: exp
    type(TpetraMap) :: tgt
    integer(size_type), allocatable :: n(:)
    integer :: my_image_id, num_images, tmp, cond
    OUT0("Starting TpetraImportExport_getNumExportIDs!")

    my_image_id = comm%getRank()
    num_images = comm%getSize()

    if (num_images < 2) return;

    ! Code here lifted from previous test
    if (my_image_id == 0 .or. my_image_id == num_images-1) then
      allocate(n(2))
    else
      allocate(n(3))
    end if

    if (my_image_id == 0) then
      n(:) = [my_image_id+1, my_image_id+2]
    else if (my_image_id == num_images - 1) then
      n(:) = [my_image_id, my_image_id+1]
    else
      n(:) = [my_image_id, my_image_id+1, my_image_id+2]
    end if

    cond = merge(1, 2, (my_image_id == 0 .or. my_image_id == num_images-1))

    ! Import test
    call TPetra_IE_MakeImport_NU(comm, 1, n, imp)
    tmp = imp%getNumExportIDs()
    TEST_ASSERT(tmp == cond)

    ! Export test
    call TPetra_IE_MakeExport_NU(comm, n, 1, exp)
    tmp = exp%getNumExportIDs()
    TEST_ASSERT(tmp == cond)

    deallocate(n)
    call imp%release(); TEST_IERR()
    call exp%release(); TEST_IERR()
    call tgt%release(); TEST_IERR()

    OUT0("Finished TpetraImportExport_getNumExportIDs!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_getNumExportIDs)

 ! ------------------------------getSourceMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_getSourceMap)
    type(TpetraImport) :: imp
    type(TpetraExport) :: exp
    type(TpetraMap) :: src, comp

    OUT0("Starting TpetraImportExport_getSourceMap!")

    ! Make both maps contiguous and uniform (1-to-1 required for import and export)
    src = TpetraMap(TPETRA_GLOBAL_INVALID, 1 * comm%getSize(), comm)

    ! Import test
    call TPetra_IE_MakeImport(comm, 1 * comm%getSize(), 5 * comm%getSize(), imp)
    comp = imp%getSourceMap()
    TEST_ASSERT(src%isSameAs(comp))

    ! Export test
    call comp%release(); TEST_IERR()
    call TPetra_IE_MakeExport(comm, 1 * comm%getSize(), 5 * comm%getSize(), exp)
    comp = exp%getSourceMap()
    TEST_ASSERT(src%isSameAs(comp))

    call comp%release(); TEST_IERR()
    call src%release(); TEST_IERR()
    call imp%release(); TEST_IERR()
    call exp%release(); TEST_IERR()

    OUT0("Finished TpetraImportExport_getSourceMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_getSourceMap)

 ! ------------------------------getTargetMap------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImportExport_getTargetMap)
    type(TpetraImport) :: imp
    type(TpetraExport) :: exp
    type(TpetraMap) :: tgt, comp

    OUT0("Starting TpetraImportExport_getTargetMap!")

    ! Make both maps contiguous and uniform (1-to-1 required for import and export)
    tgt = TpetraMap(TPETRA_GLOBAL_INVALID, 5 * comm%getSize(), comm)

    ! Import test
    call TPetra_IE_MakeImport(comm, 1 * comm%getSize(), 5 * comm%getSize(), imp)
    comp = imp%getTargetMap()
    TEST_ASSERT(tgt%isSameAs(comp))

    ! Export test
    call comp%release(); TEST_IERR()
    call TPetra_IE_MakeExport(comm, 1 * comm%getSize(), 5 * comm%getSize(), exp)
    comp = exp%getTargetMap()
    TEST_ASSERT(tgt%isSameAs(comp))

    call comp%release(); TEST_IERR()
    call tgt%release(); TEST_IERR()
    call imp%release(); TEST_IERR()
    call exp%release(); TEST_IERR()

    OUT0("Finished TpetraImportExport_getTargetMap!")

  END_FORTRILINOS_UNIT_TEST(TpetraImportExport_getTargetMap)

#if 0
  !Obsolete test
  ! ----------------------------isLocallyComplete----------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImport_isLocallyComplete)
    type(TpetraImport) :: imp
    type(TpetraMap) :: src, tgt
    integer :: exps, i
    integer, allocatable :: expPIDs(:)
    OUT0("Starting TpetraImport_isLocallyComplete!")

    src = TpetraMap(TPETRA_GLOBAL_INVALID, 1, comm)
    tgt = TpetraMap(TPETRA_GLOBAL_INVALID, 1, comm)

    imp = TpetraImport(src, tgt); TEST_IERR()
    exps = imp%getNumExportIDs(); TEST_IERR()
    allocate(expPIDs(exps)))
    expPIDs = imp%getExportPIDs(); TEST_IERR()
    do i = 1, size(expPIDs)
       TEST_ASSERT(expPIDs(i) /= -1)
    end do

    deallocate(expPIDs)
    call imp%release();   TEST_IERR()
    call src%release();   TEST_IERR()
    call tgt%release();   TEST_IERR()

    OUT0("Finished TpetraImport_isLocallyComplete!")

  END_FORTRILINOS_UNIT_TEST(TpetraImport_isLocallyComplete)
#endif

  ! --------------------------createRemoteOnlyImport-------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraImport_createRemoteOnlyImport)
    type(TpetraImport) :: importer, remote_only_importer
    type(TpetraMap) :: smap, dmap, remote_only_dmap
    integer(size_type) :: remote_id_count, expected_id_count
    integer :: gid, insert_gid_index, s_count, d_count, num_procs
    integer(global_ordinal_type), allocatable :: gids(:)

    OUT0("Starting TpetraImport_createRemoteOnlyImport!")

    num_procs = comm%getSize()

    ! only tests in parallel - note Trilinos test also skips testing
    ! createRemoteOnlyImport is num_procs = 1.
    ! it appeasrs getNumRemoteIDs() gives the wrong value for serial or MPI 1 rank
    ! it should be 0 but will include dmap elements not in smap
    if (num_procs == 1) then
      OUT0("Should not be calling TpetraImport_createRemoteOnlyImport because there is only 1 rank.")
      return
    end if

    ! Note test is setup to pass for any values here
    s_count = 5
    d_count = 10

    smap = TpetraMap(TPETRA_GLOBAL_INVALID, s_count, comm); TEST_IERR()
    dmap = TpetraMap(TPETRA_GLOBAL_INVALID, d_count, comm); TEST_IERR()

    ! For example, 4 ranke:
    ! Rank    0       1       2       3
    ! smap   1-5      6-10    11-15   16-20
    ! dmap   1-10     11-20   21-30   31-40

    ! Now we need a new map which contains only the remote entries
    ! That means entries in dmap which are in smap but not local
    !        6-10     11-20    none    none

    importer = TpetraImport(smap, dmap); TEST_IERR()
    remote_id_count = importer%getNumRemoteIDs()

    ! Make the new map with the proper count
    allocate(gids(remote_id_count))

    ! This can be eeplaced with getRemoteLIDs or used to test getRemoteLIDs
    ! We woulld set each line like this:
    !  do k = 1, remote_id_count
    !    gids(k) = dmap%getGlobalElement(importer%getRemoteLIDs()[k])
    !  end do
    ! Currently getRemoteLIDs is ignored due to +- 1 issues to resolved.

    ! Manually determine expected_id_count to compare
    ! Also fill gids
    ! We only care about ids in both, so loop 1 to min of s_count, d_count * num_procs
    ! So we say it's in the remote map if:
    !       it's in both maps
    !   and it's not local in the smap
    !   and it's local in the dmap
    expected_id_count = 0
    do gid = 1, min(d_count, s_count) * num_procs
      ! Check if gid is not local to our smap - below or above our min/max local range
      if ((gid <= s_count * comm%getRank()) .or. (gid > s_count * comm%getRank() + s_count)) then
         ! Check if it is local to our dmap
         if ((gid > d_count * comm%getRank()) .and. (gid <= d_count * comm%getRank() + d_count)) then
          gids(expected_id_count+1) = gid
          expected_id_count = expected_id_count + 1
        end if
      end if
    end do

    ! Validate the total count was filled as expected
    TEST_EQUALITY(remote_id_count, expected_id_count)

    remote_only_dmap = TpetraMap(TPETRA_GLOBAL_INVALID, gids, comm); TEST_IERR()

    remote_only_importer = importer%createRemoteOnlyImport(remote_only_dmap); TEST_IERR()

    TEST_EQUALITY(remote_id_count, remote_only_importer%getNumRemoteIDs())

    call importer%release();   TEST_IERR()
    call remote_only_importer%release();   TEST_IERR()
    call smap%release();   TEST_IERR()
    call dmap%release();   TEST_IERR()
    call remote_only_dmap%release();   TEST_IERR()
    deallocate(gids)

    OUT0("Finished TpetraImport_createRemoteOnlyImport!")

  END_FORTRILINOS_UNIT_TEST(TpetraImport_createRemoteOnlyImport)
end program test_TpetraImportExport
