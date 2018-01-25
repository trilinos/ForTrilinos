! Copyright 2017, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
program main
#include "ForTrilinosTpetra_config.hpp"

#ifdef HAVE_MPI
  use mpi
#endif
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra

  implicit none
  type(TeuchosComm) :: comm
  integer(global_size_type), parameter :: invalid = -1
  integer :: i
  type(TpetraMap) :: Map
  type(TpetraMultiVector) :: A, B, A2
  integer(size_type), parameter :: num_vecs=2, num_local=2, LDA=2
  real(scalar_type), parameter :: zero=0., one=1., two=2., four=4., negone=-1.
  real(scalar_type) :: norms(num_vecs), zeros(num_vecs), values(6)

#ifdef HAVE_MPI
  call MPI_INIT(i)
  call comm%create(MPI_COMM_WORLD)
#else
  call comm%create()
#endif

  zeros = zero

  call Map%create(invalid, num_local, comm)

  values = [one, one, two, two, four, four]
  call A%create(map, values(1:4), LDA, num_vecs)

  call A2%create(A, TeuchosCopy)
  call A2%release()

  call A2%create(A, TeuchosCopy)
  call A2%release()

  call A%release()
  call Map%release()
  call comm%release()

end program main
