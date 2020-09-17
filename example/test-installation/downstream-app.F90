! Copyright 2020-, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program downstream_app
! --------------------------------------------------------------------------- !
! Example downstream application.
! --------------------------------------------------------------------------- !
#include "ForTrilinos_config.h"
#include "ForTrilinos.h"

use iso_fortran_env
use, intrinsic :: iso_c_binding

use forerror, only : get_fortrilinos_version
use forteuchos
use fortpetra

#if FORTRILINOS_USE_MPI
use mpi
#endif

implicit none

! -- ForTrilinos objects
type(TeuchosComm) :: comm

! -- Scalars
integer :: ierr
integer :: my_rank

! --------------------------------------------------------------------------- !


! Initialize MPI subsystem, if applicable
#if FORTRILINOS_USE_MPI
call MPI_INIT(ierr)
if (ierr /= 0) then
  stop "MPI failed to init"
endif
comm = TeuchosComm(MPI_COMM_WORLD)
#else
comm = TeuchosComm()
#endif

my_rank = comm%getRank()
write(*,*) "Hello from rank", my_rank

if (my_rank == 0) then
  write(*,*) "Using ForTrilinos version: ", get_fortrilinos_version()
endif

#if FORTRILINOS_USE_MPI
! Finalize MPI must be called after releasing all handles
call MPI_FINALIZE(ierr)
#endif

end program downstream_app
