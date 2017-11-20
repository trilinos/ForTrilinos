program test_TpetraMap
#include "FortranTestMacros.h"
#include "ForTrilinosTpetra_config.hpp"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra
#ifdef HAVE_MPI
use mpi
#endif
  implicit none
  integer ierr2
#ifdef HAVE_MPI
! Initialize MPI subsystem
type(TeuchosComm) :: comm
call MPI_INIT(ierr)
EXPECT_EQ(ierr, 0)
call comm%create(MPI_COMM_WORLD)
#else
call comm%create()
#endif
  ierr2 = 0
  ierr2 = test_isOneToOne()
  if (ierr2 /= 0) then
    write(*,*) "Test 'isOneToOne' failed!"
  end if
  ierr2 = test_getGlobalNumElements()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getGlobalNumElements' failed!"
  end if
  ierr2 = test_getNodeNumElements()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getNodeNumElements' failed!"
  end if
  ierr2 = test_getIndexBase()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getIndexBase' failed!"
  end if
  ierr2 = test_getMinLocalIndex()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getMinLocalIndex' failed!"
  end if
  ierr2 = test_getMaxLocalIndex()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getMaxLocalIndex' failed!"
  end if
  ierr2 = test_getMinGlobalIndex()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getMinGlobalIndex' failed!"
  end if
  ierr2 = test_getMaxGlobalIndex()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getMaxGlobalIndex' failed!"
  end if
  ierr2 = test_getMinAllGlobalIndex()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getMinAllGlobalIndex' failed!"
  end if
  ierr2 = test_getMaxAllGlobalIndex()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getMaxAllGlobalIndex' failed!"
  end if
  ierr2 = test_getLocalElement()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getLocalElement' failed!"
  end if
  ierr2 = test_getGlobalElement()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getGlobalElement' failed!"
  end if
  ierr2 = test_getNodeElementList()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getNodeElementList' failed!"
  end if
  ierr2 = test_isNodeLocalElement()
  if (ierr2 /= 0) then
    write(*,*) "Test 'isNodeLocalElement' failed!"
  end if
  ierr2 = test_isNodeGlobalElement()
  if (ierr2 /= 0) then
    write(*,*) "Test 'isNodeGlobalElement' failed!"
  end if
  ierr2 = test_isUniform()
  if (ierr2 /= 0) then
    write(*,*) "Test 'isUniform' failed!"
  end if
  ierr2 = test_isContiguous()
  if (ierr2 /= 0) then
    write(*,*) "Test 'isContiguous' failed!"
  end if
  ierr2 = test_isDistributed()
  if (ierr2 /= 0) then
    write(*,*) "Test 'isDistributed' failed!"
  end if
  ierr2 = test_isCompatible()
  if (ierr2 /= 0) then
    write(*,*) "Test 'isCompatible' failed!"
  end if
  ierr2 = test_isSameAs()
  if (ierr2 /= 0) then
    write(*,*) "Test 'isSameAs' failed!"
  end if
  ierr2 = test_locallySameAs()
  if (ierr2 /= 0) then
    write(*,*) "Test 'locallySameAs' failed!"
  end if
  ierr2 = test_getComm()
  if (ierr2 /= 0) then
    write(*,*) "Test 'getComm' failed!"
  end if
  ierr2 = test_description()
  if (ierr2 /= 0) then
    write(*,*) "Test 'description' failed!"
  end if
  ierr2 = test_removeEmptyProcesses()
  if (ierr2 /= 0) then
    write(*,*) "Test 'removeEmptyProcesses' failed!"
  end if
  ierr2 = test_replaceCommWithSubset()
  if (ierr2 /= 0) then
    write(*,*) "Test 'replaceCommWithSubset' failed!"
  end if
  EXPECT_EQ(ierr2, 0)
  call comm%release()
#ifdef HAVE_MPI
    ! Finalize MPI must be called after releasing all handles
    call MPI_FINALIZE(ierr2)
    EXPECT_EQ(0, ierr2)
#endif
contains

! ---------------------------------isOneToOne--------------------------------- !
  integer function test_isOneToOne()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global, index_base, indices(4)
    ierr2 = 0
    index_base = 1
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm)
    if (.not. Obj%isOneToOne()) then
        write(*,*) "isOneToOne: Expected map to be one to one"
        ierr2 = 1
    end if
    call Obj%release()
    if (comm%getSize() > 1) then
        indices = [1, 2, 3, 4]
        call Obj%create(num_global, indices, index_base, comm)
        if (Obj%isOneToOne()) then
            write(*,*) "isOneToOne: Expected map to NOT be one to one"
            ierr2 = 1
        end if
        call Obj%release()
    end if
    test_isOneToOne = ierr2
    return
  end function

! ----------------------------getGlobalNumElements---------------------------- !
  integer function test_getGlobalNumElements()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    index_base = 1
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getGlobalNumElements()
    if (fresult /= num_global) then
        write(*,*) "getGlobalNumElements: Expected ", num_global, " elements, got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getGlobalNumElements = ierr2
    return
  end function

! -----------------------------getNodeNumElements----------------------------- !
  integer function test_getNodeNumElements()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_SIZE_T) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    index_base = 1
    num_global = 4*comm%getSize()
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getNodeNumElements()
    if (fresult /= 4) then
        write(*,*) "getNodeNumElements: Expected ", 4, " elements, got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getNodeNumElements = ierr2
    return
  end function

! --------------------------------getIndexBase-------------------------------- !
  integer function test_getIndexBase()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getIndexBase()
    if (fresult /= 1) then
        write(*,*) "getIndexBase: Expected indexBase = ", 1, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    index_base = 0
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getIndexBase()
    if (fresult /= 0) then
        write(*,*) "getIndexBase: Expected indexBase = ", 0, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getIndexBase = ierr2
    return
  end function

! ------------------------------getMinLocalIndex------------------------------ !
  integer function test_getMinLocalIndex()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getMinLocalIndex()
    if (fresult /= 1) then
        write(*,*) "getMinLocalIndex: Expected minLocalIndex = ", 1, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getMinLocalIndex = ierr2
    return
  end function

! ------------------------------getMaxLocalIndex------------------------------ !
  integer function test_getMaxLocalIndex()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getMaxLocalIndex()
    if (fresult /= 4) then
        write(*,*) "getMaxLocalIndex: Expected maxLocalIndex = ", 4, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getMaxLocalIndex = ierr2
    return
  end function

! -----------------------------getMinGlobalIndex------------------------------ !
  integer function test_getMinGlobalIndex()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getMinGlobalIndex()
    expected = comm%getRank() * 4 + 1
    if (fresult /= expected) then
        write(*,*) "getMindGlobalIndex: Expected minGlobalIndex = ", expected, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getMinGlobalIndex = ierr2
    return
  end function

! -----------------------------getMaxGlobalIndex------------------------------ !
  integer function test_getMaxGlobalIndex()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getMaxGlobalIndex()
    expected = comm%getRank() * 4 + 4
    if (fresult /= expected) then
        write(*,*) "getMaxGlobalIndex: Expected maxGloblIndex = ", expected, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getMaxGlobalIndex = ierr2
  end function

! ----------------------------getMinAllGlobalIndex---------------------------- !
  integer function test_getMinAllGlobalIndex()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getMinAllGlobalIndex()
    if (fresult /= index_base) then
        write(*,*) "getMinAllGlobalIndex: Expected minAllGlobalIndex = ", index_base, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getMinAllGlobalIndex = ierr2
    return
  end function

! ----------------------------getMaxAllGlobalIndex---------------------------- !
  integer function test_getMaxAllGlobalIndex()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getMaxAllGlobalIndex()
    if (fresult /= num_global) then
        write(*,*) "getMaxAllGlobalIndex: Expected maxAllGlobalIndex = ", num_global, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getMaxAllGlobalIndex = ierr2
    return
  end function

! ------------------------------getLocalElement------------------------------- !
  integer function test_getLocalElement()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    globalindex = comm%getRank() * 4 + 1
    fresult = Obj%getLocalElement(globalindex)
    if (fresult /= 1) then
        write(*,*) "Expected local element = ", 1, " got ", fresult
        ierr2 = 1
    end if
    globalindex = comm%getRank() * 4 + 4
    fresult = Obj%getLocalElement(globalindex)
    if (fresult /= 4) then
        write(*,*) "Expected local element = ", 4, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getLocalElement = ierr2
    return
  end function

! ------------------------------getGlobalElement------------------------------ !
  integer function test_getGlobalElement()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    localindex = 1
    expected = comm%getRank() * 4 + 1
    fresult = Obj%getGlobalElement(localindex)
    if (fresult /= expected) then
        write(*,*) "Expected global element = ", expected, " got ", fresult
        ierr2 = 1
    end if
    localindex = 4
    expected = comm%getRank() * 4 + 4
    fresult = Obj%getGlobalElement(localindex)
    if (fresult /= expected) then
        write(*,*) "Expected global element = ", expected, " got ", fresult
        ierr2 = 1
    end if
    call Obj%release()
    test_getGlobalElement = ierr2
    return
  end function

! -----------------------------getNodeElementList----------------------------- !
  integer function test_getNodeElementList()
    integer :: ierr2
    type(TpetraMap) :: Obj
    type(TeuchosArrayViewLongLongConst) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getNodeElementList()
    ! TODO: @aprokop: how to test fresult?  I need to check the values, but
    ! TODO: ArrayView (and friends) don't have a method to get values
    ! TODO: (that I know of)
    write(*,*) "getNodeElementList: Test not yet implemented"
    call fresult%release()
    call Obj%release()
    test_getNodeElementList = ierr2
    return
  end function

! -----------------------------isNodeLocalElement----------------------------- !
  integer function test_isNodeLocalElement()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    localindex = 1
    fresult = Obj%isNodeLocalElement(localindex)
    if (.not. fresult) then
        write(*,*) "isNodeLocalElement: expected 1 to be a local index"
        ierr2 = 1
    end if
    localindex = 5
    fresult = Obj%isNodeLocalElement(localindex)
    if (fresult) then
        write(*,*) "isNodeLocalElement: expected 5 to NOT be a local index"
        ierr2 = 1
    end if
    call Obj%release()
    test_isNodeLocalElement = ierr2
    return
  end function

! ----------------------------isNodeGlobalElement----------------------------- !
  integer function test_isNodeGlobalElement()
    integer :: ierr2
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    globalindex = int(comm%getRank() * 4 + 1, kind=global_ordinal_type)
    fresult = Obj%isNodeGlobalElement(globalindex)
    if (.not. fresult) then
        write(*,*) "isNodeGlobalElement: expected 1 to be a global index"
        ierr2 = 1
    end if
    globalindex = num_global + 1
    fresult = Obj%isNodeGlobalElement(globalindex)
    if (fresult) then
        write(*,*) "isNodeGlobalElement: expected 5 to NOT be a global index"
        ierr2 = 1
    end if
    call Obj%release()
    test_isNodeGlobalElement = ierr2
    return
  end function

! ---------------------------------isUniform---------------------------------- !
  integer function test_isUniform()
    integer :: ierr2
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer :: k
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%isUniform()
    if (.not. fresult) then
        write(*,*) "isUniform: expected map to be uniform"
        ierr2 = 1
    end if
    call Obj%release()
    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm);
    fresult = Obj%isUniform()
    if (fresult) then
        write(*,*) "isUniform: expected map to NOT be uniform"
        ierr2 = 1
    end if
    call Obj%release()
    test_isUniform = ierr2
    return
  end function

! --------------------------------isContiguous-------------------------------- !
  integer function test_isContiguous()
    integer :: ierr2
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer :: k
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%isContiguous()
    if (.not. fresult) then
        write(*,*) "isContiguous: expected map to be uniform"
        ierr2 = 1
    end if
    call Obj%release()
    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm);
    fresult = Obj%isContiguous()
    if (fresult) then
        write(*,*) "isContiguous: expected map to NOT be uniform"
        ierr2 = 1
    end if
    call Obj%release()
    test_isContiguous = ierr2
  end function

! -------------------------------isDistributed-------------------------------- !
  integer function test_isDistributed()
    integer :: ierr2
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%isDistributed()
    if (comm%getSize() == 1) then
      if (fresult) then
        write(*,*) "isDistributed: expected map to NOT be distributed"
        ierr2 = 1
      end if
    else
      if (.not. fresult) then
        write(*,*) "isDistributed: expected map to be distributed"
        ierr2 = 1
      end if
    end if
    call Obj%release()

    ! All elements have 4 entries and map has only 4 entries so the map should
    ! not be distributed
    num_global = 4
    elements = [1, 2, 3, 4]
    call Obj%create(num_global, elements, index_base, comm)
    fresult = Obj%isDistributed()
    if (comm%getSize() == 1) then
      if (fresult) then
        write(*,*) "isDistributed: expected map to NOT be distributed"
        ierr2 = 1
      end if
    else
      if (fresult) then
        write(*,*) "isDistributed: expected map to NOT be distributed"
        ierr2 = 1
      end if
    end if
    call Obj%release()
    test_isDistributed = ierr2
    return
  end function

! --------------------------------isCompatible-------------------------------- !
  integer function test_isCompatible()
    integer :: ierr2
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    integer :: k
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj1%create(num_global, index_base, comm)
    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj2%create(num_global, elements, index_base, comm);
    ! Cyclic map should be compatible
    fresult = Obj1%isCompatible(Obj2)
    if (.not. fresult) then
        ierr = 1
        write(*,*) "Expected maps to be compatible"
    end if
    call Obj2%release()
    num_global = num_global + 10
    call Obj2%create(num_global, index_base, comm);
    fresult = Obj1%isCompatible(Obj2)
    if (fresult) then
        ierr = 1
        write(*,*) "Expected maps to NOT be compatible"
    end if
    call Obj2%release()
    call Obj1%release()
    test_isCompatible = ierr2
    return
  end function

! ----------------------------------isSameAs---------------------------------- !
  integer function test_isSameAs()
    integer :: ierr2
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    logical(C_BOOL) :: fresult
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    integer(global_size_type) :: invalid
    integer :: k
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj1%create(num_global, index_base, comm)
    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj2%create(num_global, elements, index_base, comm);
    ! Cyclic map should not be SameAs
    fresult = Obj1%isSameAs(Obj2)
    if (fresult) then
        ierr = 1
        write(*,*) "Expected maps to NOT be same"
    end if
    call Obj2%release()

    invalid = -1; num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm)
    fresult = Obj1%isSameAs(Obj2)
    if (.not. fresult) then
        ierr = 1
        write(*,*) "Expected maps to be same"
    end if
    call Obj1%release()
    call Obj2%release()
    test_isSameAs = ierr2
    return
  end function

! -------------------------------locallySameAs-------------------------------- !
  integer function test_locallySameAs()
    integer :: ierr2
    type(TpetraMap) :: Obj1, Obj2
    logical(C_BOOL) :: fresult
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    integer(global_size_type) :: invalid
    integer :: k
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj1%create(num_global, index_base, comm)
    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj2%create(num_global, elements, index_base, comm);
    ! Cyclic map should not be locallySameAs
    fresult = Obj1%locallySameAs(Obj2)
    if (fresult) then
        ierr = 1
        write(*,*) "locallySameAs: Expected maps to NOT be same"
    end if
    call Obj2%release()
    invalid = -1; num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm)
    fresult = Obj1%locallySameAs(Obj2)
    if (.not. fresult) then
        ierr = 1
        write(*,*) "locallySameAs: Expected maps to be same"
    end if
    call Obj1%release()
    call Obj2%release()
    test_locallySameAs = ierr2
    return
  end function

! ----------------------------------getComm----------------------------------- !
  integer function test_getComm()
    integer :: ierr2
    type(TpetraMap) :: Obj
    type(TeuchosComm) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    ierr2 = 0
    num_global = 4*comm%getSize()
    index_base = 1
    call Obj%create(num_global, index_base, comm)
    fresult = Obj%getComm()
    ! TODO: How to test fresult, I check for the same rank and size,
    ! TODO: but that is all.  Is it enough?
    if (fresult%getRank() /= comm%getRank()) then
        ierr2 = 1
        write(*,*) "getComm: expected ranks to be same"
    end if
    if (fresult%getSize() /= comm%getSize()) then
        ierr2 = 1
        write(*,*) "getComm: expected sizes to be same"
    end if
    call Obj%release()
    test_getComm = ierr2
    return
  end function

! --------------------------------description--------------------------------- !
  integer function test_description()
    integer :: ierr2
    ierr2 = 0
    ! TODO: Implement this test?
    write(*,*) 'description: Test not yet implemented'
    test_description = ierr2
    return
  end function

! ----------------------------removeEmptyProcesses---------------------------- !
  integer function test_removeEmptyProcesses()
    integer :: ierr2
    ierr2 = 0
    ! TODO: Implement this test?
    write(*,*) 'removeEmptyProcesses: Test is not yet Implemented'
    test_removeEmptyProcesses = ierr2
    return
  end function

! ---------------------------replaceCommWithSubset---------------------------- !
  integer function test_replaceCommWithSubset()
    integer :: ierr2
    ierr2 = 0
    ! TODO: Implement this test?
    write(*,*) 'replaceCommWithSubset: Test is not yet implemented'
    test_replaceCommWithSubset = ierr2
    return
  end function

end program test_TpetraMap
