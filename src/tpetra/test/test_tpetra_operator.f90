! Copyright 2018, UT-Battelle, LLC
!
! SPDX-License-Identifier: BSD-3-Clause
! License-Filename: LICENSE
module myoperators
  use forteuchos
  use fortpetra
  implicit none

  type, extends(TpetraOperator) :: MyOperator
    contains
    procedure :: apply => apply_my_operator
  end type
contains
  subroutine apply_my_operator(self, x, y, mode, alpha, beta)
    use, intrinsic :: ISO_C_BINDING
    implicit none
    class(MyOperator), intent(in) :: self
    class(TpetraMultiVector), intent(in) :: x
    class(TpetraMultiVector), intent(inout) :: y
    integer(kind(TeuchosETransp)), intent(in) :: mode
    real(C_DOUBLE), intent(in) :: alpha
    real(C_DOUBLE), intent(in) :: beta

    ! TODO

  end subroutine
end module

program test_TpetraOperator
#include "ForTrilinosTpetra_config.hpp"
#include "FortranTestUtilities.h"
  use iso_fortran_env
  use, intrinsic :: iso_c_binding
  use forteuchos
  use fortpetra
  use myoperators

  implicit none
  type(TeuchosComm) :: comm
  integer(global_size_type), parameter :: invalid = -1
  character(len=256), parameter :: FILENAME="test_tpetra_operator.f90"


  SETUP_TEST()

#ifdef HAVE_MPI
  comm = TeuchosComm(MPI_COMM_WORLD); FORTRILINOS_CHECK_IERR()
#else
  comm = TeuchosComm(); FORTRILINOS_CHECK_IERR()
#endif

  ADD_SUBTEST_AND_RUN(TpetraOperator_Build)

  call comm%release()

  TEARDOWN_TEST()

contains

! --------------------------------description--------------------------------- !
  FORTRILINOS_UNIT_TEST(TpetraOperator_Build)
    class(TpetraOperator), allocatable :: op

    allocate(op, source=MyOperator())
    call init_TpetraOperator(op); TEST_IERR()

    call op%release()
    deallocate(op)
  END_FORTRILINOS_UNIT_TEST(TpetraOperator_Build)

end program test_TpetraOperator
