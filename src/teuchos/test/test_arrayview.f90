program main

#include "FortranTestMacros.h"
#include "ForTrilinosTeuchos_config.hpp"

  use ISO_FORTRAN_ENV
  implicit none

  call test_array()
contains

  subroutine test_array()
    use ISO_FORTRAN_ENV
    use, intrinsic :: ISO_C_BINDING
    use forteuchos
    implicit none

    integer :: sz
    integer, allocatable, dimension(:) :: farr
    type(TeuchosArrayViewInt) :: arrview_int

    sz = 10

    allocate(farr(sz))

    farr(1) = 1
    farr(2) = 2

    call arrview_int%create(farr, sz)
    EXPECT_EQ(0, ierr)

    EXPECT_EQ(arrview_int%size(), sz)
    EXPECT_EQ(0, ierr)

    ! arrview_int(0) = 2
    ! EXPECT_EQ(arrview_int(1), 1)
    ! EXPECT_EQ(arrview_int(2), 2)

    call arrview_int%release()
  end subroutine

end program

!-----------------------------------------------------------------------------!
! end of parameterlist/test.f90
!-----------------------------------------------------------------------------!
