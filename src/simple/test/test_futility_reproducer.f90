
! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program main

#include "ForTrilinosSimpleInterface_config.hpp"
#include "ForTrilinos.h"

  use ISO_FORTRAN_ENV
  use, intrinsic :: ISO_C_BINDING

#ifdef HAVE_MPI
  use mpi
#endif

#include "ForTrilinos.h"
  use fortrilinos
  use forteuchos
  use fortpetra
  use fortest

  implicit none

  integer(int_type) :: my_rank, num_procs

  integer(global_size_type) :: n_global
  integer(size_type) :: n, max_entries_per_row, lda, num_eigen = 1
  integer(int_type) :: num_eigen_int
  integer(int_type) :: row_nnz

  integer :: ierr
  integer(local_ordinal_type) :: i
  integer(global_ordinal_type) :: offset
  real(scalar_type) :: one = 1.0

  type(TeuchosComm) :: comm
  type(ParameterList) :: plist
  type(TrilinosEigenSolver) :: eigen_handle
  type(TpetraMap) :: map
  type(TpetraCrsMatrix) :: A,B
  type(TpetraMultiVector) :: X
  type(TpetraWriter) :: Writer

  real(scalar_type), dimension(:), allocatable :: evalues, evectors
  integer(global_ordinal_type), dimension(:) :: cols(2)
  real(scalar_type), dimension(:) :: vals(2)

  n = 30

#ifdef HAVE_MPI
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

  my_rank = comm%getRank()
  num_procs = comm%getSize()

  write(*,*) "Processor ", my_rank, " of ", num_procs

  ! Read in the parameterList
  plist = ParameterList("Anasazi"); FORTRILINOS_CHECK_IERR()
  call load_from_xml(plist, "davidson.xml"); FORTRILINOS_CHECK_IERR()

  num_eigen_int = num_eigen
  call plist%set("NumEV", num_eigen_int); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Step 0: Construct A and B matrices
  n_global = -1
  map = TpetraMap(n_global, n, comm); FORTRILINOS_CHECK_IERR()

  max_entries_per_row = 30
  A = TpetraCrsMatrix(map, map, max_entries_per_row, TpetraDynamicProfile)
  cols(1) = 1
  cols(2) = 2
  vals(1) = 0.0015
  vals(2) = 0.325
  call A%insertGlobalValues(INT(1,global_ordinal_type), cols(1:2), vals(1:2)); FORTRILINOS_CHECK_IERR()
  cols(1) = 2
  vals(1) = 1e-14
  call A%insertGlobalValues(INT(2,global_ordinal_type), cols(1:1), vals(2:2)); FORTRILINOS_CHECK_IERR()

  B = TpetraCrsMatrix(map, map, max_entries_per_row, TpetraDynamicProfile)
  cols(1) = 1
  vals(1) = 0.1208
  call B%insertGlobalValues(INT(1,global_ordinal_type), cols(1:1), vals(1:1)); FORTRILINOS_CHECK_IERR()
#if 0
  cols(1) = 1
  cols(2) = 2
  vals(1) = -0.117
  vals(2) = 0.184
#else
  cols(1) = 1
  cols(2) = 2
  vals(1) = -0.0
  vals(2) = 0.184
#endif
  call B%insertGlobalValues(INT(2,global_ordinal_type), cols(1:2), vals(1:2)); FORTRILINOS_CHECK_IERR()

  vals(1) = 1e-14
  vals(2) = 1.0
  do i=3,30
    cols(1)=i
    call A%insertGlobalValues(INT(i,global_ordinal_type), cols(1:1), vals(1:1)); FORTRILINOS_CHECK_IERR()
    call B%insertGlobalValues(INT(i,global_ordinal_type), cols(1:1), vals(2:2)); FORTRILINOS_CHECK_IERR()
  enddo

  call A%fillComplete(); FORTRILINOS_CHECK_IERR()
  call B%fillComplete(); FORTRILINOS_CHECK_IERR()

  Writer = TpetraWriter()
  call Writer%writeSparseFile('A.mm', A)
  call Writer%writeSparseFile('B.mm', B)

  ! The solution
  allocate(evalues(num_eigen))
  X = TpetraMultiVector(map, num_eigen)

  ! Step 0: crate a handle
  eigen_handle = TrilinosEigenSolver(); FORTRILINOS_CHECK_IERR()

  ! ------------------------------------------------------------------
  ! Explicit setup and solve
  ! ------------------------------------------------------------------
  ! Step 1: initialize a handle
  call eigen_handle%init(comm); FORTRILINOS_CHECK_IERR()

  ! Step 2: setup the problem
  call eigen_handle%setup_matrix(A); FORTRILINOS_CHECK_IERR()
  call eigen_handle%setup_matrix_rhs(B); FORTRILINOS_CHECK_IERR()
  ! Step 3: setup the solver
  call eigen_handle%setup_solver(plist); FORTRILINOS_CHECK_IERR()

  ! Step 4: solve the system
  call eigen_handle%solve(evalues, X); FORTRILINOS_CHECK_IERR()

  ! FIXME: Check the solution
  if (size(evalues) /= num_eigen) then
    write(error_unit, '(A)') 'The number of returned eigenvalues does not match!'
    stop 1
  end if
  write(*,*) "Computed eigenvalues: ", evalues(1)

  ! Step 5: clean up
  call eigen_handle%finalize(); FORTRILINOS_CHECK_IERR()

  call eigen_handle%release(); FORTRILINOS_CHECK_IERR()
  call plist%release(); FORTRILINOS_CHECK_IERR()
  call X%release(); FORTRILINOS_CHECK_IERR()
  call A%release(); FORTRILINOS_CHECK_IERR()
  call map%release(); FORTRILINOS_CHECK_IERR()
  call comm%release(); FORTRILINOS_CHECK_IERR()
  deallocate(evalues)

#ifdef HAVE_MPI
  ! Finalize MPI must be called after releasing all handles
  call MPI_FINALIZE(ierr)
#endif


end program
