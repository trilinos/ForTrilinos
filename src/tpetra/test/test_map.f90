program main

#include "FortranTestMacros.h"
#include "ForTrilinosTpetra_config.hpp"

  use ISO_FORTRAN_ENV
  implicit none

  call test_comm()
contains

  subroutine test_comm()
    use ISO_FORTRAN_ENV
    use, intrinsic :: ISO_C_BINDING
    use forteuchos
    use fortpetra
#ifdef HAVE_MPI
    use mpi
#endif
    implicit none

    type(TeuchosComm) :: comm
    type(TpetraMap) :: map
    integer(C_LONG) :: num_global, sum_local
    integer(C_SIZE_T) :: num_local

#ifdef HAVE_MPI
    ! Initialize MPI subsystem
    call MPI_INIT(ierr)
    if (ierr /= 0) then
      write(*,*) "MPI failed to init"
      stop 1
    endif

    call comm%create(MPI_COMM_WORLD)
#else
    call comm%create()
#endif

    num_local = 10 + comm%getRank()
    sum_local = 10*comm%getSize() + (comm%getSize()*(comm%getSize()-1))/2

    ! Test 1
    num_global = sum_local

    call map%create(num_global, 0, comm)
    EXPECT_EQ(ierr, 0)

    EXPECT_EQ(map%getGlobalNumElements(), num_global)

    call map%release()

    ! Test 2
    num_global = -1

    call map%create(num_global, num_local, 0, comm)
    EXPECT_EQ(ierr, 0)

    EXPECT_EQ(map%getNodeNumElements(), num_local)
    EXPECT_EQ(map%getGlobalNumElements(), sum_local)

    call map%release()


    call comm%release()
#ifdef HAVE_MPI
    ! Finalize MPI must be called after releasing all handles
    call MPI_FINALIZE(ierr)
    EXPECT_EQ(0, ierr)
#endif

  end subroutine

end program

!-----------------------------------------------------------------------------!
! end of parameterlist/test.f90
!-----------------------------------------------------------------------------!
