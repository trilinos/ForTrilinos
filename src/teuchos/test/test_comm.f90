program main

#include "FortranTestMacros.h"
#include "ForTrilinosTeuchos_config.hpp"

  use ISO_FORTRAN_ENV
  implicit none

  call test_comm()
contains

  subroutine test_comm()
    use ISO_FORTRAN_ENV
    use, intrinsic :: ISO_C_BINDING
    use forteuchos
#ifdef HAVE_MPI
    use mpi
#endif
    implicit none

    type(TeuchosComm) :: comm
    integer :: comm_rank_f, comm_rank_c
    integer :: comm_size_f, comm_size_c

#ifdef HAVE_MPI
    ! Initialize MPI subsystem
    call MPI_INIT(ierr)
    if (ierr /= 0) then
      write(*,*) "MPI failed to init"
      stop 1
    endif

    call comm%create(MPI_COMM_WORLD)
    EXPECT_EQ(0, ierr)
    EXPECT_TRUE(c_associated(comm%swigptr))

    call MPI_COMM_RANK(MPI_COMM_WORLD, comm_rank_f, ierr)
    call MPI_COMM_SIZE(MPI_COMM_WORLD, comm_size_f, ierr)
#else
    call comm%create()
    EXPECT_EQ(0, ierr)
    EXPECT_TRUE(c_associated(comm%swigptr))

    comm_rank_f = 0
    comm_size_f = 1
#endif

    comm_rank_c = comm%getRank()
    EXPECT_EQ(0, ierr)
    EXPECT_EQ(comm_rank_f, comm_rank_c)

    comm_size_c = comm%getSize()
    EXPECT_EQ(0, ierr)
    EXPECT_EQ(comm_size_f, comm_size_c)

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
