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
  integer kerr
  type(TeuchosComm) :: comm

#ifdef HAVE_MPI
  ! Initialize MPI subsystem
  call MPI_INIT(ierr)
  EXPECT_EQ(ierr, 0)
  call comm%create(MPI_COMM_WORLD)
#else
  call comm%create()
#endif

  kerr = 0

  if (test_isOneToOne() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'isOneToOne' FAILED!"
  end if

  if (test_getGlobalNumElements() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getGlobalNumElements' FAILED!"
  end if

  if (test_getNodeNumElements() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getNodeNumElements' FAILED!"
  end if

  if (test_getIndexBase() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getIndexBase' FAILED!"
  end if

  if (test_getMinLocalIndex() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getMinLocalIndex' FAILED!"
  end if

  if (test_getMaxLocalIndex() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getMaxLocalIndex' FAILED!"
  end if

  if (test_getMinGlobalIndex() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getMinGlobalIndex' FAILED!"
  end if

  if (test_getMaxGlobalIndex() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getMaxGlobalIndex' FAILED!"
  end if

  if (test_getMinAllGlobalIndex() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getMinAllGlobalIndex' FAILED!"
  end if

  if (test_getMaxAllGlobalIndex() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getMaxAllGlobalIndex' FAILED!"
  end if

  if (test_getLocalElement() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getLocalElement' FAILED!"
  end if

  if (test_getGlobalElement() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getGlobalElement' FAILED!"
  end if

  if (test_getNodeElementList() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getNodeElementList' FAILED!"
  end if

  if (test_isNodeLocalElement() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'isNodeLocalElement' FAILED!"
  end if

  if (test_isNodeGlobalElement() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'isNodeGlobalElement' FAILED!"
  end if

  if (test_isUniform() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'isUniform' FAILED!"
  end if

  if (test_isContiguous() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'isContiguous' FAILED!"
  end if

  if (test_isDistributed() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'isDistributed' FAILED!"
  end if

  if (test_isCompatible() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'isCompatible' FAILED!"
  end if

  if (test_isSameAs() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'isSameAs' FAILED!"
  end if

  if (test_locallySameAs() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'locallySameAs' FAILED!"
  end if

  if (test_getComm() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'getComm' FAILED!"
  end if

  if (test_description() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'description' FAILED!"
  end if

  if (test_removeEmptyProcesses() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'removeEmptyProcesses' FAILED!"
  end if

  if (test_replaceCommWithSubset() /= 0) then
    kerr = kerr + 1
    if (comm%getRank() == 0) &
      write(*,*) "Test 'replaceCommWithSubset' FAILED!"
  end if

  if (comm%getRank() == 0) then
    if (kerr == 0) then
      write(*,*) "Test PASSED"
    else
      write(*,*) "A total of ", kerr, " tests FAILED"
    end if
  end if

  EXPECT_EQ(kerr, 0)

  call comm%release()
#ifdef HAVE_MPI
    ! Finalize MPI must be called after releasing all handles
    call MPI_FINALIZE(ierr)
    EXPECT_EQ(0, ierr)
#endif

contains

! ---------------------------------isOneToOne--------------------------------- !
  integer function test_isOneToOne()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(global_ordinal_type) :: num_global, index_base, indices(4)
    jerr = 0
    index_base = 1
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_isOneToOne = jerr + ierr
      ierr = 0
      return
    end if

    if (.not. Obj%isOneToOne()) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isOneToOne: Expected map to be one to one"
    end if

    call Obj%release()

    if (comm%getSize() > 1) then
        indices = [1, 2, 3, 4]
        call Obj%create(num_global, indices, index_base, comm)
        if (ierr /= 0) then
          test_isOneToOne = jerr + ierr
          ierr = 0
          return
        end if

        if (Obj%isOneToOne()) then
          jerr = jerr + 1
          if (comm%getRank() == 0) &
            write(*,*) "isOneToOne: Expected map to NOT be one to one"
        end if

        call Obj%release()

    end if

    test_isOneToOne = jerr + ierr
    ierr = 0
    return

  end function

! ----------------------------getGlobalNumElements---------------------------- !
  integer function test_getGlobalNumElements()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    index_base = 1
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getGlobalNumElements = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getGlobalNumElements()
    if (fresult /= num_global) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getGlobalNumElements: Expected ", num_global, " elements, got ", fresult
    end if

    call Obj%release()

    test_getGlobalNumElements = jerr + ierr
    ierr = 0
    return

  end function

! -----------------------------getNodeNumElements----------------------------- !
  integer function test_getNodeNumElements()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_SIZE_T) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    index_base = 1
    num_global = 4*comm%getSize()

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getNodeNumElements = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getNodeNumElements()
    if (fresult /= 4) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getNodeNumElements: Expected ", 4, " elements, got ", fresult
    end if

    call Obj%release()

    test_getNodeNumElements = jerr + ierr
    ierr = 0
    return

  end function

! --------------------------------getIndexBase-------------------------------- !
  integer function test_getIndexBase()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getIndexBase = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getIndexBase()
    if (fresult /= 1) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getIndexBase: Expected indexBase = ", 1, " got ", fresult
    end if

    call Obj%release()

    index_base = 0
    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getIndexBase = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getIndexBase()
    if (fresult /= 0) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getIndexBase: Expected indexBase = ", 0, " got ", fresult
    end if

    call Obj%release()

    test_getIndexBase = jerr + ierr
    ierr = 0
    return

  end function

! ------------------------------getMinLocalIndex------------------------------ !
  integer function test_getMinLocalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getMinLocalIndex = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getMinLocalIndex()
    if (fresult /= 1) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMinLocalIndex: Expected minLocalIndex = ", 1, " got ", fresult
    end if

    call Obj%release()

    test_getMinLocalIndex = jerr + ierr
    ierr = 0
    return

  end function

! ------------------------------getMaxLocalIndex------------------------------ !
  integer function test_getMaxLocalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getMaxLocalIndex = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getMaxLocalIndex()
    if (fresult /= 4) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMaxLocalIndex: Expected maxLocalIndex = ", 4, " got ", fresult
    end if

    call Obj%release()

    test_getMaxLocalIndex = jerr + ierr
    ierr = 0
    return

  end function

! -----------------------------getMinGlobalIndex------------------------------ !
  integer function test_getMinGlobalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getMinGlobalIndex = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getMinGlobalIndex()
    expected = comm%getRank() * 4 + 1
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMindGlobalIndex: Expected minGlobalIndex = ", expected, " got ", fresult
    end if
    call Obj%release()
    test_getMinGlobalIndex = jerr + ierr
    ierr = 0
    return

  end function

! -----------------------------getMaxGlobalIndex------------------------------ !
  integer function test_getMaxGlobalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getMaxGlobalIndex = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getMaxGlobalIndex()
    expected = comm%getRank() * 4 + 4
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMaxGlobalIndex: Expected maxGloblIndex = ", expected, " got ", fresult
    end if

    call Obj%release()

    test_getMaxGlobalIndex = jerr + ierr
    ierr = 0
    return

  end function

! ----------------------------getMinAllGlobalIndex---------------------------- !
  integer function test_getMinAllGlobalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getMinAllGlobalIndex = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getMinAllGlobalIndex()
    if (fresult /= index_base) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMinAllGlobalIndex: Expected minAllGlobalIndex = ", index_base, " got ", fresult
    end if

    call Obj%release()

    test_getMinAllGlobalIndex = jerr + ierr
    ierr = 0
    return

  end function

! ----------------------------getMaxAllGlobalIndex---------------------------- !
  integer function test_getMaxAllGlobalIndex()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getMaxAllGlobalIndex = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getMaxAllGlobalIndex()
    if (fresult /= num_global) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getMaxAllGlobalIndex: Expected maxAllGlobalIndex = ", num_global, " got ", fresult
    end if

    call Obj%release()

    test_getMaxAllGlobalIndex = jerr + ierr
    ierr = 0
    return

  end function

! ------------------------------getLocalElement------------------------------- !
  integer function test_getLocalElement()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    integer(C_INT) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getLocalElement = jerr + ierr
      ierr = 0
      return
    end if

    globalindex = comm%getRank() * 4 + 1
    fresult = Obj%getLocalElement(globalindex)
    if (fresult /= 1) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getLocalElement: Expected local element = ", 1, " got ", fresult
    end if

    globalindex = comm%getRank() * 4 + 4
    fresult = Obj%getLocalElement(globalindex)
    if (fresult /= 4) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getLocalElement: Expected local element = ", 4, " got ", fresult
    end if

    call Obj%release()

    test_getLocalElement = jerr + ierr
    ierr = 0
    return

  end function

! ------------------------------getGlobalElement------------------------------ !
  integer function test_getGlobalElement()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    integer(C_LONG_LONG) :: fresult, expected
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getGlobalElement = jerr + ierr
      ierr = 0
      return
    end if

    localindex = 1
    expected = comm%getRank() * 4 + 1
    fresult = Obj%getGlobalElement(localindex)
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getGlobalElement: Expected global element = ", expected, " got ", fresult
    end if

    localindex = 4
    expected = comm%getRank() * 4 + 4
    fresult = Obj%getGlobalElement(localindex)
    if (fresult /= expected) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "getGlobalElement: Expected global element = ", expected, " got ", fresult
    end if

    call Obj%release()

    test_getGlobalElement = jerr + ierr
    ierr = 0
    return

  end function

! -----------------------------getNodeElementList----------------------------- !
  integer function test_getNodeElementList()
    integer :: jerr
    type(TpetraMap) :: Obj
    type(TeuchosArrayViewLongLongConst) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getNodeElementList = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%getNodeElementList()
    if (ierr /= 0) then
      test_getNodeElementList = jerr + ierr
      ierr = 0
      return
    end if

    ! TODO: getNodeElementList should be modified to be a subroutine that takes
    ! TODO: the element list as a return argument.  Otherwise, there is no
    ! TODO: (current) way to get to the data in the ArrayView returned.
    if (comm%getRank() == 0) &
      write(*,*) "getNodeElementList: Test not yet implemented"

    call fresult%release()
    call Obj%release()
    test_getNodeElementList = jerr + ierr
    ierr = 0
    return

  end function

! -----------------------------isNodeLocalElement----------------------------- !
  integer function test_isNodeLocalElement()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_INT) :: localindex
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_isNodeLocalElement = jerr + ierr
      ierr = 0
      return
    end if

    localindex = 1
    fresult = Obj%isNodeLocalElement(localindex)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeLocalElement: expected 1 to be a local index"
    end if

    localindex = 5
    fresult = Obj%isNodeLocalElement(localindex)
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeLocalElement: expected 5 to NOT be a local index"
    end if

    call Obj%release()

    test_isNodeLocalElement = jerr + ierr
    ierr = 0
    return

  end function

! ----------------------------isNodeGlobalElement----------------------------- !
  integer function test_isNodeGlobalElement()
    integer :: jerr
    type(TpetraMap) :: Obj
    integer(C_LONG_LONG) :: globalindex
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_isNodeGlobalElement = jerr + ierr
      ierr = 0
      return
    end if

    globalindex = int(comm%getRank() * 4 + 1, kind=global_ordinal_type)
    fresult = Obj%isNodeGlobalElement(globalindex)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeGlobalElement: expected 1 to be a global index"
    end if
    globalindex = num_global + 1

    fresult = Obj%isNodeGlobalElement(globalindex)
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isNodeGlobalElement: expected 5 to NOT be a global index"
    end if

    call Obj%release()

    test_isNodeGlobalElement = jerr + ierr
    ierr = 0
    return

  end function

! ---------------------------------isUniform---------------------------------- !
  integer function test_isUniform()
    integer :: jerr
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer :: k
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_isUniform = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%isUniform()
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isUniform: expected map to be uniform"
    end if

    call Obj%release()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm);
    if (ierr /= 0) then
      test_isUniform = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%isUniform()
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isUniform: expected map to NOT be uniform"
    end if

    call Obj%release()

    test_isUniform = jerr + ierr
    ierr = 0
    return

  end function

! --------------------------------isContiguous-------------------------------- !
  integer function test_isContiguous()
    integer :: jerr
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer :: k
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_isContiguous = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%isContiguous()
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isContiguous: expected map to be uniform"
    end if

    call Obj%release()

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj%create(num_global, elements, index_base, comm);
    if (ierr /= 0) then
      test_isContiguous = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%isContiguous()
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank() == 0) &
        write(*,*) "isContiguous: expected map to NOT be uniform"
    end if
    call Obj%release()
    test_isContiguous = jerr + ierr
  end function

! -------------------------------isDistributed-------------------------------- !
  integer function test_isDistributed()
    integer :: jerr
    type(TpetraMap) :: Obj
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_isDistributed = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%isDistributed()
    if (comm%getSize() == 1) then

      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isDistributed: expected map to NOT be distributed"
      end if

    else

      if (.not. fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isDistributed: expected map to be distributed"
      end if

    end if

    call Obj%release()

    ! All elements have 4 entries and map has only 4 entries so the map should
    ! not be distributed
    num_global = 4
    elements = [1, 2, 3, 4]
    call Obj%create(num_global, elements, index_base, comm)
    if (ierr /= 0) then
      test_isDistributed = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj%isDistributed()
    if (comm%getSize() == 1) then

      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isDistributed: expected map to NOT be distributed"
      end if

    else

      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isDistributed: expected map to NOT be distributed"
      end if

    end if

    call Obj%release()

    test_isDistributed = jerr + ierr
    ierr = 0
    return

  end function

! --------------------------------isCompatible-------------------------------- !
  integer function test_isCompatible()
    integer :: jerr
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    logical(C_BOOL) :: fresult
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    integer :: k
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj1%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_isCompatible = jerr + ierr
      ierr = 0
      return
    end if

    do k = 1, 4
      elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
    end do
    call Obj2%create(num_global, elements, index_base, comm);
    if (ierr /= 0) then
      test_isCompatible = jerr + ierr
      ierr = 0
      return
    end if

    ! Cyclic map should be compatible
    fresult = Obj1%isCompatible(Obj2)
    if (.not. fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isCompatible: Expected maps to be compatible"
    end if

    call Obj2%release()

    num_global = num_global + 10
    call Obj2%create(num_global, index_base, comm);
    if (ierr /= 0) then
      test_isCompatible = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj1%isCompatible(Obj2)
    if (fresult) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "isCompatible: Expected maps to NOT be compatible"
    end if

    call Obj2%release()
    call Obj1%release()

    test_isCompatible = jerr + ierr
    ierr = 0
    return

  end function

! ----------------------------------isSameAs---------------------------------- !
  integer function test_isSameAs()
    integer :: jerr
    type(TpetraMap) :: Obj1, Obj2
    type(TpetraMap) :: map
    logical(C_BOOL) :: fresult
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    integer(global_size_type) :: invalid
    integer :: k
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj1%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_isSameAs = jerr + ierr
      ierr = 0
      return
    end if

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      call Obj2%create(num_global, elements, index_base, comm);
      if (ierr /= 0) then
        test_isSameAs = jerr + ierr
        ierr = 0
        return
      end if

      ! Cyclic map should not be SameAs
      fresult = Obj1%isSameAs(Obj2)
      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "isSameAs: Expected maps to NOT be same"
      end if

      call Obj2%release()

    end if

    invalid = -1; num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm)
    if (ierr /= 0) then
      test_isSameAs = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj1%isSameAs(Obj2)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "isSameAs: Expected maps to be same"
    end if

    call Obj1%release()
    call Obj2%release()

    test_isSameAs = jerr + ierr
    ierr = 0
    return

  end function

! -------------------------------locallySameAs-------------------------------- !
  integer function test_locallySameAs()
    integer :: jerr
    type(TpetraMap) :: Obj1, Obj2
    logical(C_BOOL) :: fresult
    integer(size_type) :: num_local
    integer(global_ordinal_type) :: num_global, index_base, elements(4)
    integer(global_size_type) :: invalid
    integer :: k
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj1%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_locallySameAs = jerr + ierr
      ierr = 0
      return
    end if

    if (comm%getSize() > 1) then

      do k = 1, 4
        elements(k) = int(comm%getRank()+k*comm%getSize(), kind=global_ordinal_type)
      end do
      call Obj2%create(num_global, elements, index_base, comm);
      if (ierr /= 0) then
        test_locallySameAs = jerr + ierr
        ierr = 0
        return
      end if

      ! Cyclic map should not be locallySameAs
      fresult = Obj1%locallySameAs(Obj2)
      if (fresult) then
        jerr = jerr + 1
        if (comm%getRank()==0) &
          write(*,*) "locallySameAs: Expected maps to NOT be same"
      end if

      call Obj2%release()

    end if

    invalid = -1; num_local = 4
    call Obj2%create(invalid, num_local, index_base, comm)
    if (ierr /= 0) then
      test_locallySameAs = jerr + ierr
      ierr = 0
      return
    end if

    fresult = Obj1%locallySameAs(Obj2)
    if (.not. fresult) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "locallySameAs: Expected maps to be same"
    end if

    call Obj1%release()
    call Obj2%release()

    test_locallySameAs = jerr + ierr
    ierr = 0
    return

  end function

! ----------------------------------getComm----------------------------------- !
  integer function test_getComm()
    integer :: jerr
    type(TpetraMap) :: Obj
    type(TeuchosComm) :: fresult
    integer(global_ordinal_type) :: num_global, index_base
    jerr = 0
    num_global = 4*comm%getSize()
    index_base = 1

    call Obj%create(num_global, index_base, comm)
    if (ierr /= 0) then
      test_getComm = jerr + ierr
      ierr = 0
      return
    end if

    ! We only test the comm returned has the same rank and size.  More
    ! comprehensive testing is (assumed to be) done in Teuchos itself.
    fresult = Obj%getComm()
    if (fresult%getRank() /= comm%getRank()) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "getComm: expected ranks to be same"
    end if

    if (fresult%getSize() /= comm%getSize()) then
      jerr = jerr + 1
      if (comm%getRank()==0) &
        write(*,*) "getComm: expected sizes to be same"
    end if

    call Obj%release()

    test_getComm = jerr + ierr
    ierr = 0
    return

  end function

! --------------------------------description--------------------------------- !
  integer function test_description()
    integer :: jerr
    jerr = 0
    ! TODO: Implement this test?
    if (comm%getRank()==0) &
      write(*,*) 'description: Test not yet implemented'
    test_description = jerr + ierr
    ierr = 0
    return
  end function

! ----------------------------removeEmptyProcesses---------------------------- !
  integer function test_removeEmptyProcesses()
    integer :: jerr
    jerr = 0
    ! TODO: Implement this test?
    if (comm%getRank()==0) &
      write(*,*) 'removeEmptyProcesses: Test is not yet Implemented'
    test_removeEmptyProcesses = jerr + ierr
    ierr = 0
    return
  end function

! ---------------------------replaceCommWithSubset---------------------------- !
  integer function test_replaceCommWithSubset()
    integer :: jerr
    jerr = 0
    ! TODO: Implement this test?
    if (comm%getRank()==0) &
      write(*,*) 'replaceCommWithSubset: Test is not yet implemented'
    test_replaceCommWithSubset = jerr + ierr
    ierr = 0
    return
  end function

end program test_TpetraMap
