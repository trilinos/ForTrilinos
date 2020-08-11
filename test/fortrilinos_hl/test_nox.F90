! Copyright 2017-2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
! ---------------------------------------------------------------------------- !

program main
  ! -------------------------------------------------------------------------- !
#include "ForTrilinos_config.h"
#include "ForTrilinos.h"
  use ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING
  use fortrilinos_hl
  use forteuchos
#if FORTRILINOS_USE_MPI
  use mpi
#endif
  implicit none
  integer :: ierr
  type(TeuchosComm) :: comm
  type(ParameterList) :: params
  ! -------------------------------------------------------------------------- !

#if FORTRILINOS_USE_MPI
  ! Initialize MPI subsystem
  call MPI_INIT(ierr)
  if (ierr /= 0) then
    write(*,*) "MPI failed to init"
    stop 1
  endif
  comm = TeuchosComm(MPI_COMM_WORLD)
#else
  comm = TeuchosComm()
#endif

  params = ParameterList("TpetraModelEvaluator1DFEM")
  call load_from_xml(params, 'nox_params.xml')
  call main2(comm, params, ierr)

  call params%release(); FORTRILINOS_CHECK_IERR()
  call comm%release(); FORTRILINOS_CHECK_IERR()

  if (ierr /= 0) then
    write(*,*) "End Result: TEST FAILED"
    stop 1
  else
    write(*,*) "End Result: TEST PASSED"
  endif

contains

  ! -------------------------------------------------------------------------- !

  subroutine main2(comm, params, ierr)
    ! ------------------------------------------------------------------------ !
    use fortpetra
    use fortrilinos_hl
    use TpetraModelEvaluator1DFEM_module
    implicit none
    integer :: ierr
    type(TeuchosComm) :: comm
    type(ParameterList) :: params
    class(ForModelEvaluator), allocatable :: evaluator
    type(NOXSolver) :: nox_solver
    integer(global_size_type) :: num_global_elems
    real(scalar_type) :: z_min, z_max
    integer(kind(NOXStatusType)) :: status
    ! ------------------------------------------------------------------------ !

    ierr = 0

    ! Create the model evaluator object
    num_global_elems = 100
    z_min = 0.0
    z_max = 1.0
    allocate(evaluator, source=TpetraModelEvaluator1DFEM(comm, num_global_elems, z_min, z_max))

    call init_ForModelEvaluator(evaluator); FORTRILINOS_CHECK_IERR()
    call evaluator%setup(params); FORTRILINOS_CHECK_IERR()

    nox_solver = NOXSolver(evaluator)
    call nox_solver%setup(params); FORTRILINOS_CHECK_IERR()
    status = nox_solver%solve(); FORTRILINOS_CHECK_IERR()

    if (status /= NOXConverged) ierr = 1

    call evaluator%release(); FORTRILINOS_CHECK_IERR()
    deallocate(evaluator)

    call params%release(); FORTRILINOS_CHECK_IERR()
    call nox_solver%release(); FORTRILINOS_CHECK_IERR()

  end subroutine main2

end program
